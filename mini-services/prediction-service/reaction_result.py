"""ReactionResult - 从论文源码改编的反应结果处理类

去除了Django依赖（prediction_run对象），改为接受纯RDKit参数。
保留所有核心功能：
- 严格反应中心识别（get_strict_reaction_centres）
- 2D构象优化（minimise_products_2d_conformation）
- 3D结构能量最小化（minimise_products_3d_structure）
- 立体化学检测（detect_sterochemistry_from_3d_structure）
- 最大键长计算（get_max_bond_length）
- 键变化检测（get_bond_changes）
"""

import logging
from typing import Collection, Dict, List, Optional, Set, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdMolTransforms
from rdkit.Geometry.rdGeometry import Point2D, Point3D

logger = logging.getLogger("ezmechanism")


class ReactionResult:
    """用于存储和操作化学反应结果的类（论文源码改编版）

    本类实现了EzMechanism论文中描述的关键功能：
    1. 严格反应中心识别 - 通过比较反应前后键的对称差集精确识别参与反应的原子
    2. 3D结构优化 - 使用MMFF94s力场对产物进行能量最小化
    3. 立体化学检测 - 从3D构象中检测和更新立体化学信息
    4. 键长评估 - 计算新形成键的长度用于筛选不合理产物

    与论文原始版本的区别：
    - 移除了对Django prediction_run对象的依赖
    - 所有方法使用纯RDKit参数
    """

    # 日志提示信息
    LOG_DISCARD_BOND_CARBON_RESIDUES = "与催化残基的碳原子形成了键"
    LOG_DISCARD_UNSTABLE_PRODUCTS = "产物包含不稳定的子结构"
    LOG_DISCARD_SANITIZE_ERROR = "无法 Sanitize 该产物："
    LOG_DISCARD_REACTION_WITH_R_GROUP = "该产物涉及 R 基团"
    LOG_DISCARD_ONLY_INERT_ATOMS_OF_RXN_SUBSTRATE = "该产物仅涉及底物中的惰性原子"

    def __init__(
        self,
        reactants: Collection[Chem.Mol],
        products: Collection[Chem.Mol],
        rule_smarts: str = "",
        am_to_3d_coord: Optional[Dict[int, 'Point3D']] = None,
        new_bonds_from_rule: Optional[set] = None,
    ):
        """
        初始化 ReactionResult 对象

        Upload upload_reaction_result.py L23-43:
        self.reactants = reactants
        self.rule = rule
        self.predict_mechanism = prediction_run
        self.products, self.max_bond_length, self.reaction_centres_am = self._quick_parse(products)

        参数：
        :param reactants: 反应物分子集合（需要有atom map编号）
        :param products: RDKit RunReactants生成的产物分子列表
        :param rule_smarts: 反应规则SMARTS字符串（用于标识规则）
        :param am_to_3d_coord: 反应物PDB 3D坐标映射（Upload: self.predict_mechanism.am_to_3d_coord）
        :param new_bonds_from_rule: 规则定义的新键集合（Upload: self.predict_mechanism.rule_to_bond_changes_am_heavy[self.rule.id][0]）
        """
        self.reactants = list(reactants)
        self.products = list(products)
        self.rule_smarts = rule_smarts

        # 反应中心原子映射编号集合（由get_strict_reaction_centres填充）
        self.reaction_centres_am: Set[int] = set()

        # 新形成的键列表（由get_bond_changes填充）
        # 格式: [(frozenset((am1, am2)), bond_type), ...]
        self.new_bonds: List[tuple] = []

        # 断裂的键列表
        # 格式: [(frozenset((am1, am2)), bond_type), ...]
        self.broken_bonds: List[tuple] = []

        # Upload 格式的新键列表（由get_bond_changes填充，供get_max_bond_length使用）
        # 格式: [(product, p_idx1, p_idx2, reactant1, r_idx1, reactant2, r_idx2), ...]
        self.new_bonds_detailed: List[tuple] = []

        # 最大键长
        # _quick_parse 设置为 PDB 快速估算值（Dijkstra 评分用）
        # get_max_bond_length（步骤⑦）覆盖为精确值（画图用）
        self.max_bond_length: float = 0.0

        # 规则映射编号到产物映射编号的映射关系
        self.rule_am_to_product_am: Dict[int, int] = {}

        # Upload 外部数据（对齐 Upload: self.predict_mechanism.am_to_3d_coord / rule_to_bond_changes_am_heavy）
        self.am_to_3d_coord = am_to_3d_coord or {}
        self._new_bonds_from_rule = new_bonds_from_rule

        # Upload: self.products, self.max_bond_length, self.reaction_centres_am = self._quick_parse(products)
        self._quick_parse(products)

    def _quick_parse(self, products):
        """
        快速解析产物，提取用于优先级计算的关键信息

        Upload upload_reaction_result.py L45-84 (完整复制):
          1. 将产物中的 old_mapno 属性转换为原子映射编号
          2. 更新原子映射编号（用 isotope 保存的 reactant AM 替换 rule template AM）
          3. 基于 PDB 3D 坐标计算新形成键的最大长度（Dijkstra 评分用）
          4. 提取反应中心的原子映射编号

        Upload upload_reaction_result.py L58-84:
          new_bonds_from_rule = self.predict_mechanism.rule_to_bond_changes_am_heavy[self.rule.id][0]
          ...
          for bond in new_bonds_from_rule:
              atom1_am, atom2_am = bond
              bond_length = (self.predict_mechanism.am_to_3d_coord[self.rule_am_to_product_am[atom1_am]] -
                             self.predict_mechanism.am_to_3d_coord[self.rule_am_to_product_am[atom2_am]]).Length()
              max_bond_length = max(max_bond_length, bond_length)
          return products, max_bond_length, rcs_ams

        Zero-tolerance coordinate policy (aligned with Upload source code):
          In PDB mode (am_to_3d_coord is provided and non-empty), if any atom in
          the new bonds lacks a PDB coordinate, raise KeyError to discard this
          reaction result entirely. This matches the Upload source code behavior
          where am_to_3d_coord always contains all atoms — a missing key means
          the intermediate doesn't have real PDB coordinates and the result is
          invalid for Dijkstra scoring.

        返回值：
        :return: (products, max_bond_length, rcs_ams)
                  max_bond_length 为 PDB 快速估算值，用于 Dijkstra 评分
        :raises KeyError: In PDB mode, when an atom in new bonds lacks PDB coordinates
        """
        # Upload L58-59: new_bonds_from_rule = self.predict_mechanism.rule_to_bond_changes_am_heavy[self.rule.id][0]
        new_bonds_from_rule = self._new_bonds_from_rule
        self.rule_am_to_product_am = {}
        rcs_ams: Set[int] = set()

        for product in products:
            for atom in product.GetAtoms():
                if atom.GetAtomicNum() > 1 and atom.HasProp("old_mapno"):
                    rule_am = atom.GetIntProp("old_mapno")
                    atom.ClearProp("old_mapno")
                    product_am = atom.GetIsotope()
                    self.rule_am_to_product_am[rule_am] = product_am
                    atom.SetAtomMapNum(product_am)
                    rcs_ams.add(product_am)

        # Upload L75-82: 基于 PDB 坐标计算新形成键的最大长度（Dijkstra 评分用）
        #
        # Upload 系统中 am_to_3d_coord 始终有值（来自 PDB 结构），
        # 所有原子的坐标在初始化时就已全部进入字典。
        # 源代码直接用 self.am_to_3d_coord[am] 访问，缺失则 KeyError → 丢弃整个结果。
        #
        # 我们的对齐策略：PDB 模式下，缺失坐标 → 抛出 KeyError → 丢弃该反应结果。
        # 这样保证了评分仅基于真实 PDB 键长，杜绝伪坐标污染。
        max_bond_length = 0
        if new_bonds_from_rule and self.am_to_3d_coord:
            # PDB mode: zero tolerance — missing coordinates → KeyError → discard
            for bond in new_bonds_from_rule:
                atom1_am, atom2_am = bond
                p_am1 = self.rule_am_to_product_am.get(atom1_am)
                p_am2 = self.rule_am_to_product_am.get(atom2_am)
                if p_am1 is None or p_am2 is None:
                    continue
                # Source code style: direct dict access → KeyError if missing
                bond_length = (self.am_to_3d_coord[p_am1] -
                               self.am_to_3d_coord[p_am2]).Length()
                max_bond_length = max(max_bond_length, bond_length)

        elif new_bonds_from_rule and not self.am_to_3d_coord:
            # Non-PDB mode: no 3D coordinates available, scoring falls back to
            # the penalty mechanism in _get_reaction_data_for_products
            pass

        # Upload L84: return products, max_bond_length, rcs_ams
        self.max_bond_length = max_bond_length
        self.reaction_centres_am = rcs_ams
        self._initial_rcs_ams = rcs_ams

    def _generate_and_backfill_coords(
        self,
        products: Collection[Chem.Mol],
        missing_ams: Set[int],
    ) -> None:
        """Generate 3D conformers for product molecules that lack coordinates
        and backfill am_to_3d_coord with the generated coordinates.

        ⚠️ DEPRECATED — This method is NO LONGER called by _quick_parse.
        In PDB mode (am_to_3d_coord provided), _quick_parse now raises KeyError
        when atoms lack coordinates, which discards the entire reaction result
        (aligned with Upload source code behavior). This method is kept for
        potential future use in non-PDB mode scenarios.

        Original documentation preserved for reference:

        This is called by _quick_parse when intermediate molecules (produced by
        RunReactants) have atoms whose coordinates are not in am_to_3d_coord.
        The original Upload system doesn't need this because am_to_3d_coord
        always contains all required atoms (populated from PDB at initialization).

        ⚠️ SIDE EFFECT — Mutates self.am_to_3d_coord (and by reference,
        run_state.am_to_3d_coord):
          This method directly writes generated coordinates into the shared
          am_to_3d_coord dict. Because self.am_to_3d_coord is a reference to
          run_state.am_to_3d_coord (not a copy), all backfilled entries
          propagate to the global search state immediately. This means:
            • Positive: Subsequent Dijkstra iterations can score bonds involving
              these atoms using the newly available coordinates.
            • Risk: Generated coordinates are approximate (ETKDGv3 + MMFF, not
              PDB-derived), and errors may accumulate over multiple steps.
              However, the "only add missing" policy (PDB coords are never
              overwritten) ensures that accurate PDB coordinates always take
              priority, limiting error propagation.

        Strategy:
          1. For each product molecule that contains missing atoms:
             a. Try EmbedMolecule + MMFFOptimizeMolecule to generate a 3D conformer
             b. Extract coordinates for ALL non-H atoms with AM != 0
             c. Backfill into am_to_3d_coord (Python dict reference → propagates
                to run_state.am_to_3d_coord, so subsequent steps benefit)
          2. If EmbedMolecule fails, try ETKDGv3 as fallback
          3. PDB coordinates are never overwritten — only missing entries are added

        Args:
            products: Product molecules from RunReactants
            missing_ams: Set of atom map numbers that need coordinates
        """
        if not missing_ams or not self.am_to_3d_coord:
            # If am_to_3d_coord is empty, there's no PDB data at all —
            # generating coords won't help because Dijkstra has no reference frame.
            # The penalty score (score=9) in _get_reaction_data_for_products
            # will handle this case.
            return

        remaining = set(missing_ams)

        for product in products:
            if not remaining:
                break  # All missing coords filled

            # Check if this product contains any of the missing atoms
            product_ams = set()
            for atom in product.GetAtoms():
                if atom.GetAtomicNum() > 1:  # Non-H only
                    am = atom.GetAtomMapNum()
                    if am in remaining:
                        product_ams.add(am)

            if not product_ams:
                continue  # This product doesn't contain missing atoms

            # Try to generate a 3D conformer for this product
            try:
                mol = Chem.Mol(product)  # Deep copy to avoid modifying original
                mol = Chem.AddHs(mol)

                # Strategy 1: ETKDGv3 (best quality, may fail for some molecules)
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.numThreads = 1
                embed_result = AllChem.EmbedMolecule(mol, params)

                if embed_result == -1:
                    # Strategy 2: Basic embedding with random coordinates
                    embed_result = AllChem.EmbedMolecule(
                        mol, randomSeed=42, useRandomCoords=True
                    )

                if embed_result == -1:
                    logger.debug(
                        "EmbedMolecule failed for product (SMILES: %s), "
                        "skipping coordinate generation",
                        Chem.MolToSmiles(product, canonical=True)[:40],
                    )
                    continue

                # MMFF optimization for reasonable geometry
                try:
                    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
                except Exception:
                    pass  # Non-critical — unoptimized coords are still better than none

                # Extract coordinates from conformer 0 and backfill am_to_3d_coord
                # We match atoms by AtomMapNum (set by _quick_parse to product AM)
                conformer = mol.GetConformer()
                n_backfilled = 0
                for atom in mol.GetAtoms():
                    am = atom.GetAtomMapNum()
                    if am in remaining and atom.GetAtomicNum() > 1:
                        pos = conformer.GetAtomPosition(atom.GetIdx())
                        # Only add if not already present (PDB coords take priority)
                        if am not in self.am_to_3d_coord:
                            self.am_to_3d_coord[am] = Point3D(pos.x, pos.y, pos.z)
                            remaining.discard(am)
                            n_backfilled += 1

                if n_backfilled > 0:
                    logger.debug(
                        "_generate_and_backfill_coords: added %d atoms to "
                        "am_to_3d_coord (remaining: %d)",
                        n_backfilled, len(remaining),
                    )

            except Exception as e:
                logger.debug(
                    "3D coordinate generation failed for product: %s", str(e)[:80]
                )
                continue

        if remaining:
            logger.debug(
                "_generate_and_backfill_coords: %d atoms could not be filled",
                len(remaining),
            )

    # =========================================================================
    # 核心功能1: 严格反应中心识别
    # =========================================================================

    def get_strict_reaction_centres(self) -> Set[int]:
        """
        确定严格反应中心，即反应中发生键变化的原子

        通过比较反应物和产物中的键集合（对称差集）来精确识别
        哪些原子参与了键的形成或断裂。

        Upload upload_reaction_result.py L108-137:
        原代码仅返回集合，不覆写 self.reaction_centres_am。
        self.reaction_centres_am 由 _quick_parse 设置为规则匹配的宽集合（rcs_ams），
        供 get_bonds_away_from_rxn_rc 使用。如果覆写为严格集合，
        会导致 RC discount 和 3D 约束范围出错。

        返回值：
        :return: 严格反应中心的原子映射编号集合（仅键变化原子，不覆写 self）
        """
        reactant_bonds_ams: Set[tuple] = set()
        product_bonds_ams: Set[tuple] = set()

        # 提取反应物中的键信息（对齐 Upload：不过滤 AM=0 的键）
        for reactant in self.reactants:
            for bond in reactant.GetBonds():
                bond_atoms = (bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())
                bond_type = bond.GetBondTypeAsDouble()
                reactant_bonds_ams.add((frozenset(bond_atoms), bond_type))

        # 提取产物中的键信息（对齐 Upload：不过滤 AM=0 的键）
        for product in self.products:
            for bond in product.GetBonds():
                bond_atoms = (bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())
                bond_type = bond.GetBondTypeAsDouble()
                product_bonds_ams.add((frozenset(bond_atoms), bond_type))

        # 找出反应前后键的变化（对称差集）
        strict_reaction_centres: Set[int] = set()
        for bond_diff in reactant_bonds_ams ^ product_bonds_ams:
            strict_reaction_centres.update(bond_diff[0])

        # 对齐 Upload: 仅返回，不覆写 self.reaction_centres_am
        # self.reaction_centres_am 保持 _quick_parse 设置的宽集合（rcs_ams）
        return strict_reaction_centres

    def get_bonds_away_from_rxn_rc(
        self,
        overall_rc_atom_maps: Optional[Set[int]] = None,
    ) -> Optional[int]:
        """计算产物与整体反应中心之间的距离（键的数量）。

        对齐 Upload 的 upload_reaction_result.py L224-235:
            if self.reaction_centres_am & self.predict_mechanism.prediction.step_rc_ams:
                return 1

        在 Upload 中，step_rc_ams 来自 Django DB（管理员标注的反应中心原子映射号），
        通过 rxn_am_to_scheme_am 映射到 scheme 级别。

        在我们系统中，overall_rc_atom_maps 来自
        identify_reaction_center_atom_maps()，基于 Ketcher RXN 的 atom-mapped SMILES
        计算，等效于 paper 的 step_rc_ams。

        Args:
            overall_rc_atom_maps: 整体反应的反应中心原子映射号集合。
                即 paper 的 step_rc_ams。如果为 None 或空集，返回 None。

        Returns:
            若产物涉及反应中心（原子映射号有交集），返回 1；否则返回 None。
            （与 Upload 实现一致：暂未实现具体距离计算）
        """
        if not overall_rc_atom_maps:
            return None
        if self.reaction_centres_am & overall_rc_atom_maps:
            return 1
        return None

    # =========================================================================
    # 核心功能2: 键变化检测
    # =========================================================================

    def get_bond_changes(self) -> Tuple[List[tuple], List[tuple]]:
        """
        获取反应中形成和断裂的键信息

        同时填充 self.new_bonds（AM格式）和 self.new_bonds_detailed（Upload格式）。
        Upload格式供 get_max_bond_length 使用，包含产物/反应物分子和原子索引。

        返回值：
        :return: 两个列表，分别存储新形成的键和断裂的键
                AM格式: (frozenset((am1, am2)), bond_type)
        """
        self.new_bonds = []
        self.broken_bonds = []
        self.new_bonds_detailed = []

        # 收集所有产物原子的映射信息
        product_atom_maps: Dict[int, Tuple[Chem.Mol, int]] = {}
        for product in self.products:
            for atom in product.GetAtoms():
                am = atom.GetAtomMapNum()
                if am != 0:
                    product_atom_maps[am] = (product, atom.GetIdx())

        # 收集所有反应物原子的映射信息
        reactant_atom_maps: Dict[int, Tuple[Chem.Mol, int]] = {}
        for reactant in self.reactants:
            for atom in reactant.GetAtoms():
                am = atom.GetAtomMapNum()
                if am != 0:
                    reactant_atom_maps[am] = (reactant, atom.GetIdx())

        # 找出产物中存在但反应物中不存在的键（新形成的键）
        product_bond_set: Set[tuple] = set()
        for product in self.products:
            for bond in product.GetBonds():
                am1 = bond.GetBeginAtom().GetAtomMapNum()
                am2 = bond.GetEndAtom().GetAtomMapNum()
                if am1 != 0 and am2 != 0:
                    product_bond_set.add((frozenset((am1, am2)), bond.GetBondTypeAsDouble()))

        reactant_bond_set: Set[tuple] = set()
        for reactant in self.reactants:
            for bond in reactant.GetBonds():
                am1 = bond.GetBeginAtom().GetAtomMapNum()
                am2 = bond.GetEndAtom().GetAtomMapNum()
                if am1 != 0 and am2 != 0:
                    reactant_bond_set.add((frozenset((am1, am2)), bond.GetBondTypeAsDouble()))

        # 新形成的键 — 同时生成 Upload 格式的详细数据
        for bond_info in product_bond_set - reactant_bond_set:
            bond_ams, bond_type = bond_info
            self.new_bonds.append((bond_ams, bond_type))

            # 生成 Upload 详细格式: (product, p_idx1, p_idx2, reactant1, r_idx1, reactant2, r_idx2)
            ams = list(bond_ams)
            if len(ams) == 2:
                am1, am2 = ams
                p_info1 = product_atom_maps.get(am1)
                p_info2 = product_atom_maps.get(am2)
                r_info1 = reactant_atom_maps.get(am1)
                r_info2 = reactant_atom_maps.get(am2)
                if p_info1 and p_info2 and r_info1 and r_info2:
                    self.new_bonds_detailed.append(
                        (p_info1[0], p_info1[1], p_info2[1],
                         r_info1[0], r_info1[1], r_info2[0], r_info2[1])
                    )

        # 断裂的键
        for bond_info in reactant_bond_set - product_bond_set:
            bond_ams, bond_type = bond_info
            self.broken_bonds.append((bond_ams, bond_type))

        return self.new_bonds, self.broken_bonds

    # =========================================================================
    # 核心功能3: 2D构象优化
    # =========================================================================

    def minimise_products_2d_conformation(self):
        """
        优化产物的2D几何结构，使其更合理

        基于现有坐标构建新的2D构象，保持重原子的相对位置不变，
        重新计算最优的2D布局。
        """
        for product in self.products:
            try:
                product.UpdatePropertyCache()
                conformer_2d = product.GetConformer(0)

                # 对齐 Upload: Point2D 直接从 Point3D 构造
                coord_map = {atom.GetIdx(): Chem.rdGeometry.Point2D(conformer_2d.GetAtomPosition(atom.GetIdx()))
                             for atom in product.GetAtoms() if atom.GetSymbol() != "H"}
                new_conformer_id = Chem.rdDepictor.Compute2DCoords(
                    product, clearConfs=False, coordMap=coord_map
                )
                new_conformer = product.GetConformer(new_conformer_id)

                # 更新原子位置并移除多余构象
                for atom in product.GetAtoms():
                    conformer_2d.SetAtomPosition(
                        atom.GetIdx(), new_conformer.GetAtomPosition(atom.GetIdx())
                    )
                product.RemoveConformer(new_conformer_id)
            except Exception as e:
                logger.warning(f"2D构象优化失败: {e}")

    # =========================================================================
    # 核心功能4: 3D结构能量最小化
    # =========================================================================

    def minimise_products_3d_structure(self):
        """
        对产物的3D结构进行能量最小化，使其更稳定

        注意：
        1. 使用 MMFF94s 力场算法
        2. 对卤素原子（如 Cl）会跳过优化，因其可能导致能量最小化失败
        3. 非反应中心的重原子施加较强约束（权重5）
        4. 反应中心的重原子施加较弱约束（权重1），允许较大位移
        """
        for product in self.products:
            # 需要 confId=1（Upload 原版硬编码方式，3D构象）
            if product.GetNumConformers() < 2:
                continue

            product_copy = Chem.Mol(product)
            try:
                Chem.SanitizeMol(product_copy)
            except Exception as e:
                logger.warning(f"无法清理分子结构: {e}")
                continue

            # 获取分子属性和力场对象
            try:
                mp = Chem.AllChem.MMFFGetMoleculeProperties(
                    product_copy, mmffVariant='MMFF94s'
                )
                if mp is None:
                    continue

                ff = Chem.AllChem.MMFFGetMoleculeForceField(
                    product_copy, mp, confId=1
                )
            except Exception:
                continue

            if ff is None:
                logger.warning(f"力场无法识别分子")
                continue

            # 检查分子中是否包含卤素原子，若包含则跳过优化（对齐 Upload：仅 Cl）
            skip = False
            for atom in product_copy.GetAtoms():
                if atom.GetSymbol() in ["Cl"]:
                    logger.warning(f"跳过优化（含卤素原子）")
                    skip = True
                    break

            # 对齐 Upload: 约束在 skip 检查之前施加
            for atom in product_copy.GetAtoms():
                if atom.GetSymbol() != "H":
                    if atom.GetAtomMapNum() not in self.reaction_centres_am:
                        # 非反应中心的重原子施加较强约束
                        ff.MMFFAddPositionConstraint(atom.GetIdx(), 0.1, 5)
                    else:
                        # 反应中心的重原子施加较弱约束
                        ff.MMFFAddPositionConstraint(atom.GetIdx(), 0.1, 1)

            if not skip:
                ff.Minimize(maxIts=2000)

                # 更新产物的3D构象（confId=1）
                for atom in product.GetAtoms():
                    product.GetConformer(1).SetAtomPosition(
                        atom.GetIdx(),
                        product_copy.GetConformer(1).GetAtomPosition(atom.GetIdx()),
                    )

    # =========================================================================
    # 核心功能5: 立体化学检测
    # =========================================================================

    def detect_stereochemistry_from_3d_structure(self):
        """
        根据产物的3D结构检测并更新立体化学信息

        执行以下操作：
        1. 根据3D结构分配原子手性
        2. 清理键的手性标记（避免错误标记）
        3. 检测键的立体化学
        """
        for product in self.products:
            # 对齐 Upload: 根据 3D 结构分配原子手性
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(product, confId=1)

            # 清理键的手性标记（避免错误标记）
            for bond in product.GetBonds():
                if bond.GetStereo() == Chem.rdchem.BondStereo.STEREOANY:
                    bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)

            # 检测键的立体化学
            Chem.rdmolops.DetectBondStereoChemistry(product, product.GetConformer(1))

    # =========================================================================
    # 核心功能6: 最大键长计算（对齐 Upload upload_reaction_result.py L237-274）
    # =========================================================================

    def get_max_bond_length(self, allow_rotations: bool = True) -> float:
        """计算反应结果中的最大键长（Upload 原版方法）

        Upload upload_reaction_result.py L237-274:
        遍历 self.new_bonds（格式：product, idx1, idx2, r1, r1_idx, r2, r2_idx），
        使用反应物构象中的原子坐标计算键长，对端基原子（仅有1个键的原子）应用旋转修正。

        我们用 self.new_bonds_detailed（由 get_bond_changes 生成，格式与 Upload 的 self.new_bonds 相同）替代。
        """
        if not self.new_bonds_detailed:
            # 如果还没有计算键变化，先计算
            self.get_bond_changes()

        max_bond_length = 0.0

        for product, idx1, idx2, r1, r1_idx, r2, r2_idx in self.new_bonds_detailed:
            # 对齐 Upload: 直接使用 r_mols[i].GetConformer(1)
            correction = 0.0
            r_mols = (r1, r2)
            r_idxs = (r1_idx, r2_idx)
            pos_to_compare = []

            for i, rc_idx in enumerate((idx1, idx2)):
                r_atom = r_mols[i].GetAtomWithIdx(r_idxs[i])
                r_conformer = r_mols[i].GetConformer(1)

                if len(r_atom.GetBonds()) == 1 and allow_rotations:
                    # 端基原子：应用旋转修正
                    other_atom_idx = r_atom.GetBonds()[0].GetOtherAtom(r_atom).GetIdx()
                    bond_length = rdMolTransforms.GetBondLength(r_conformer, r_idxs[i], other_atom_idx)
                    correction -= bond_length * 0.5
                    pos_to_compare.append(r_conformer.GetAtomPosition(other_atom_idx))
                else:
                    pos_to_compare.append(r_conformer.GetAtomPosition(r_idxs[i]))

            bond_length = (pos_to_compare[0] - pos_to_compare[1]).Length() + correction
            max_bond_length = max(max_bond_length, bond_length)

        self.max_bond_length = max_bond_length
        return max_bond_length

    def compute_max_bond_length_from_3d(self, allow_rotations: bool = True) -> float:
        """向后兼容别名 — 委托给 get_max_bond_length"""
        return self.get_max_bond_length(allow_rotations=allow_rotations)

    # =========================================================================
    # 辅助功能：产物过滤
    # =========================================================================

    def has_reasonable_bond_length(self, max_bond_length: float = 6.0) -> bool:
        """
        检查产物中的新形成键是否在合理范围内

        注意：此方法当前不再用于过滤（对齐 Upload），
        键长仅作为评分参数，不做硬截止。

        参数：
        :param max_bond_length: 最大允许键长（Angstroms）
        :return: True如果所有新键都在合理范围内
        """
        if self.max_bond_length == 0:
            self.get_max_bond_length()
        return self.max_bond_length <= max_bond_length

    def has_reaction_with_r_group(self) -> bool:
        """
        检查产物是否涉及R基团（未定义的子结构）

        :return: True如果产物涉及R基团
        """
        for product in self.products:
            for atom in product.GetAtoms():
                if atom.GetSymbol() in ("R", "*") or atom.GetAtomicNum() == 0:
                    return True
        return False

    def has_bond_to_carbon_residue(self, residue_elements: Optional[Set[str]] = None) -> bool:
        """
        检查产物是否与催化残基的碳原子形成了键

        :param residue_elements: 残基元素集合
        :return: True如果产物与残基碳形成了键
        """
        if residue_elements is None:
            residue_elements = {"C", "N", "O", "S"}

        for product in self.products:
            for atom in product.GetAtoms():
                if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() == 0:
                    # 没有映射编号的重原子可能是残基原子
                    return True
        return False

    def sanitize_products(self) -> bool:
        """
        对产物进行sanitize检查

        :return: True如果所有产物都成功sanitize
        """
        for product in self.products:
            try:
                Chem.SanitizeMol(product)
            except Exception:
                return False
        return True

    # =========================================================================
    # 完整流水线：执行所有优化步骤
    # =========================================================================

    def process_full_pipeline(
        self,
        optimize_2d: bool = True,
        optimize_3d: bool = True,
        detect_stereo: bool = True,
        max_bond_length: float = 6.0,
    ) -> Dict[str, any]:
        """执行完整的产物处理流水线

        论文设计意图（破译 Upload 死代码函数的调用逻辑）：

        Upload 的搜索流程极简：
          1. ReactionResult.__init__ → _quick_parse（设置 rcs_ams）
          2. add_and_get_equivalent_mols（去重 + 价态异常检查）
          3. score 计算

        但 Upload 中定义了 5 个未被调用的方法，它们构成了完整的产物处理流水线。
        合理的调用顺序是：

          Step 1: sanitize_products()
            - 验证产物化学合理性，sanitize 失败的产物应丢弃

          Step 2: get_strict_reaction_centres()
            - 精确识别键变化原子（对称差集），用于后续步骤的约束范围

          Step 3: get_bond_changes()
            - 识别新键和断裂键，填充 new_bonds / new_bonds_detailed
            - new_bonds_detailed 供 get_max_bond_length 使用

          Step 4: minimise_products_2d_conformation()
            - 2D坐标优化，为3D优化做准备

          Step 5: minimise_products_3d_structure()
            - MMFF94s 能量最小化，使用 strict RC 原子施加弱约束

          Step 6: detect_stereochemistry_from_3d_structure()
            - 从优化后的3D构象检测立体化学

          Step 7: get_max_bond_length()
            - 基于反应物构象坐标计算最大键长（含旋转修正）→ 覆盖 self.max_bond_length（画图用）
            - 注: _quick_parse 已设置 self.max_bond_length 为 PDB 快速估算值（Dijkstra 评分用）

        当前启用的过滤（对齐 Upload 精神）：
          - Sanitize 失败 → 丢弃（对齐 Upload 的 AtomValenceException 检查）
          - 残基碳成键 → 丢弃

        已注释掉的过滤（论文设计了LOG常量但未实际实现）：
          - R基团检查
          - 键长硬截止

        参数：
        :param optimize_2d: 是否执行2D优化
        :param optimize_3d: 是否执行3D优化
        :param detect_stereo: 是否检测立体化学
        :param max_bond_length: 键长参考值（仅用于记录，不再做硬截止过滤）

        返回值：
        :return: 包含处理结果和统计信息的字典
        """
        result = {
            "success": True,
            "reaction_centres": set(),
            "new_bonds_count": 0,
            "broken_bonds_count": 0,
            "max_bond_length": 0.0,
            "bond_length_reasonable": True,
            "has_r_group": False,
            "has_residue_bond": False,
            "sanitize_ok": True,
        }

        # Step 1: Sanitize（唯一启用的过滤之一）
        if not self.sanitize_products():
            result["sanitize_ok"] = False
            result["success"] = False
            return result

        # Step 2: 严格反应中心识别（仅返回，不覆写 self.reaction_centres_am）
        result["reaction_centres"] = self.get_strict_reaction_centres()

        # Step 3: 键变化检测（填充 new_bonds + new_bonds_detailed）
        new_bonds, broken_bonds = self.get_bond_changes()
        result["new_bonds_count"] = len(new_bonds)
        result["broken_bonds_count"] = len(broken_bonds)

        # Step 3b: 残基碳成键检查（唯一启用的过滤之二）
        result["has_residue_bond"] = self.has_bond_to_carbon_residue()

        # [对齐 Upload] R-group 检查 — 注释掉，不过滤
        # 论文设计了 LOG_DISCARD_REACTION_WITH_R_GROUP，但 Upload 未实现过滤
        # result["has_r_group"] = self.has_reaction_with_r_group()

        # Step 4: 2D构象优化
        if optimize_2d:
            try:
                self.minimise_products_2d_conformation()
            except Exception:
                pass

        # Step 5: 3D结构能量最小化
        if optimize_3d:
            try:
                self.minimise_products_3d_structure()
            except Exception:
                pass

        # Step 6: 立体化学检测
        if detect_stereo:
            self.detect_stereochemistry_from_3d_structure()

        # Step 7: 最大键长计算（Upload 方法，使用反应物构象坐标）
        self.get_max_bond_length()
        result["max_bond_length"] = round(self.max_bond_length, 3)
        # [对齐 Upload] 键长仅用于评分，不做硬截止过滤
        result["bond_length_reasonable"] = True  # 始终为 True

        return result


# =========================================================================
# 独立辅助函数（从论文 predict_mechanism_util.py 改编）
# =========================================================================

def mol_is_proton(mol: Chem.Mol) -> bool:
    """检查分子是否是质子（H+）"""
    return mol.GetNumAtoms() == 1 and mol.GetAtoms()[0].GetSymbol() == "H"


def mol_is_water_species(mol: Chem.Mol) -> bool:
    """检查分子是否是水分子、氢氧根或水合氢离子"""
    return (mol.GetNumHeavyAtoms() == 1 and
            any([atom.GetSymbol() == "O" for atom in mol.GetAtoms()]))


def mol_is_hydronium(mol: Chem.Mol) -> bool:
    """检查分子是否是水合氢离子（H3O+）"""
    return mol.GetNumHeavyAtoms() == 1 and any(
        [atom.GetSymbol() == "O" and atom.GetFormalCharge() == 1
         for atom in mol.GetAtoms()]
    )


def mol_is_water_or_proton(mol: Chem.Mol) -> bool:
    """检查分子是否是水分子或质子"""
    return mol_is_proton(mol) or mol_is_water_species(mol)


def compare_mols(mol_a: Chem.Mol, mol_b: Chem.Mol,
                 atom_compare=Chem.rdFMCS.AtomCompare.CompareElements,
                 bond_compare=Chem.rdFMCS.BondCompare.CompareOrder) -> float:
    """计算两个分子之间相似原子的百分比"""
    mcs = Chem.rdFMCS.FindMCS(
        [mol_a, mol_b], atomCompare=atom_compare,
        bondCompare=bond_compare, timeout=2
    )
    if mcs.numAtoms == 0:
        return 0.0
    return mcs.numAtoms / max(mol_a.GetNumAtoms(), mol_b.GetNumAtoms())


def get_atom_map_to_point_coord(mol: Chem.Mol, conformer_id: int = 0) -> Dict[Chem.Atom, Point3D]:
    """
    return a dict that maps the atom maps of this mol to 3d coordinates

    :param mol: a rkdit Mol object with one or more conformers
    """
    try:
        conformer = mol.GetConformer(conformer_id)
    except ValueError:
        logger.warning("{} has no conformer {}".format(Chem.MolToSmiles(mol), conformer_id))
        return {}
    Point = Point3D if conformer.Is3D() else Point2D
    no_coords = 3 if conformer.Is3D() else 2
    return {atom.GetAtomMapNum(): Point(*list(conformer.GetAtomPosition(atom.GetIdx()))[:no_coords])
            for atom in mol.GetAtoms()}


def filter_products(
    products: List[Chem.Mol],
    reactants: List[Chem.Mol],
    max_bond_length: float = 6.0,
    remove_water_proton: bool = True,
) -> List[Chem.Mol]:
    """
    对产物进行多层过滤，返回合理的产物列表

    过滤规则（来自论文）：
    1. 产物无法sanitize → 丢弃
    2. 产物包含R基团 → 丢弃
    3. 产物仅涉及惰性原子 → 丢弃
    4. 新形成键过长 → 丢弃
    5. （可选）移除水和质子等小分子

    参数：
    :param products: 产物分子列表
    :param reactants: 反应物分子列表
    :param max_bond_length: 最大允许键长
    :param remove_water_proton: 是否移除水和质子

    返回值：
    :return: 过滤后的产物列表
    """
    filtered = []

    for product_tuple in products:
        valid_products = []

        for p in product_tuple:
            # Filter 1: Sanitize
            try:
                Chem.SanitizeMol(p)
            except Exception:
                logger.debug(ReactionResult.LOG_DISCARD_SANITIZE_ERROR + Chem.MolToSmiles(p))
                continue

            # Filter 2: R基团检查
            has_r = False
            for atom in p.GetAtoms():
                if atom.GetSymbol() in ("R", "*") or atom.GetAtomicNum() == 0:
                    has_r = True
                    break
            if has_r:
                logger.debug(ReactionResult.LOG_DISCARD_REACTION_WITH_R_GROUP)
                continue

            # Filter 3: 催化残基碳原子成键检查
            # 检查产物是否与催化残基的碳原子形成了键（论文中的过滤原因之一）
            has_residue_bond = False
            for atom in p.GetAtoms():
                if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() == 0:
                    # 没有映射编号的重原子可能是残基原子
                    has_residue_bond = True
                    break
            if has_residue_bond:
                logger.debug(ReactionResult.LOG_DISCARD_BOND_CARBON_RESIDUES)
                continue

            valid_products.append(p)

        if not valid_products:
            continue

        # 创建ReactionResult并执行完整流水线
        rr = ReactionResult(reactants, valid_products)
        pipeline_result = rr.process_full_pipeline(
            optimize_2d=True,
            optimize_3d=False,  # 3D优化较慢，仅在需要时启用
            detect_stereo=True,
            max_bond_length=max_bond_length,
        )

        if not pipeline_result["success"]:
            continue

        # Filter: 键长检查
        if not pipeline_result["bond_length_reasonable"]:
            continue

        # Filter: 移除水和质子
        if remove_water_proton:
            final_products = []
            for p in rr.products:
                if not mol_is_water_or_proton(p):
                    final_products.append(p)
            if not final_products:
                # 如果过滤后没有产物了，保留原始产物
                final_products = rr.products
        else:
            final_products = rr.products

        filtered.extend(final_products)

    return filtered


# =========================================================================
# 独立 Mol 优化函数（从 ReactionResult 方法提取，用于结果展示阶段）
# =========================================================================

def optimize_mol_2d(mol: Chem.Mol) -> bool:
    """对单个 Mol 做 2D 坐标优化（提取自 ReactionResult.minimise_products_2d_conformation）

    仅在 Mol 恰好只有 1 个 conformer 时调用，避免误伤 3D 构象。
    搜索产生的中间体分子通常只有 Conformer(0)（2D），优化后 SVG 渲染更美观。

    Args:
        mol: 需要优化的 RDKit Mol 对象（必须至少有 1 个 conformer）

    Returns:
        True 如果优化成功，False 如果跳过或失败
    """
    # 安全检查：仅对恰好只有 1 个 conformer 的 Mol 调用
    # NumConformers >= 2 意味着有 3D 构象 (confId=1)，不应修改 confId=0
    # NumConformers == 0 意味着无构象，无法优化
    if mol.GetNumConformers() != 1:
        return False

    try:
        mol.UpdatePropertyCache()
        conformer_2d = mol.GetConformer(0)

        # 对齐 ReactionResult: Point2D 直接从 Point3D 构造
        coord_map = {
            atom.GetIdx(): Point2D(conformer_2d.GetAtomPosition(atom.GetIdx()))
            for atom in mol.GetAtoms() if atom.GetSymbol() != "H"
        }
        new_conformer_id = Chem.rdDepictor.Compute2DCoords(
            mol, clearConfs=False, coordMap=coord_map
        )
        new_conformer = mol.GetConformer(new_conformer_id)

        # 更新原子位置并移除多余构象
        for atom in mol.GetAtoms():
            conformer_2d.SetAtomPosition(
                atom.GetIdx(), new_conformer.GetAtomPosition(atom.GetIdx())
            )
        mol.RemoveConformer(new_conformer_id)
        return True
    except Exception as e:
        logger.debug("optimize_mol_2d failed: %s", e)
        return False


def detect_mol_stereo(mol: Chem.Mol, conf_id: int = 1) -> bool:
    """对单个 Mol 做立体化学检测（提取自 ReactionResult.detect_stereochemistry_from_3d_structure）

    仅在 Mol 具有真正的 3D 构象时调用。
    搜索产生的中间体通常没有 3D 构象，此函数会安全跳过并返回 False。

    Args:
        mol: 需要检测的 RDKit Mol 对象
        conf_id: 3D 构象的 ID（默认 1，与 ReactionResult 一致）

    Returns:
        True 如果检测成功执行，False 如果跳过（无 3D 构象）或失败
    """
    # 检查构象是否存在
    if mol.GetNumConformers() <= conf_id:
        return False

    conf = mol.GetConformer(conf_id)

    # 检查构象是否真的是 3D
    # 2D 构象调用 AssignAtomChiralTagsFromStructure 会产生无意义的结果
    if not conf.Is3D():
        return False

    try:
        # 对齐 ReactionResult: 根据 3D 结构分配原子手性
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(mol, confId=conf_id)

        # 清理键的手性标记（避免错误标记）
        for bond in mol.GetBonds():
            if bond.GetStereo() == Chem.rdchem.BondStereo.STEREOANY:
                bond.SetStereo(Chem.rdchem.BondStereo.STEREONONE)

        # 检测键的立体化学
        Chem.rdmolops.DetectBondStereoChemistry(mol, conf)
        return True
    except Exception as e:
        logger.debug("detect_mol_stereo failed: %s", e)
        return False
