# EZmechanism 项目代码包 — 部署说明

> 打包时间: 2025-04-28

---

## 压缩包清单

| 编号 | 文件名 | 内容 | 大小 |
|------|--------|------|------|
| 01 | `01-prediction-service.zip` | Python 后端预测服务代码（~40个.py文件） | 核心文件 |
| 02 | `02-frontend-src.zip` | Next.js 前端源码（src/ 全部） | UI + API 路由 |
| 03 | `03-ketcher-config.zip` | Ketcher 化学编辑器配置（无node_modules） | 需 npm install |
| 04 | `04-project-config.zip` | 项目根配置（package.json, tsconfig, Caddyfile等） | 必需 |
| 05 | `05-data-database.zip` | 规则数据(xlsx) + SQLite数据库(custom.db) | 51,637条规则 |
| 06 | `06-public-assets.zip` | 静态资源（Ketcher WASM, 3Dmol.js, logo等） | 静态文件 |

---

## 本地部署步骤

### 1. 创建项目目录

```bash
mkdir ezmechanism && cd ezmechanism
```

### 2. 解压所有压缩包

```bash
# 按顺序解压（覆盖即可）
unzip 01-prediction-service.zip
unzip 02-frontend-src.zip
unzip 03-ketcher-config.zip
unzip 04-project-config.zip
unzip 05-data-database.zip
unzip 06-public-assets.zip
```

### 3. 安装前端依赖

```bash
npm install    # 或 bun install
```

### 4. 安装 Ketcher 依赖并构建

```bash
cd mini-services/ketcher-service
npm install
npm run build
# 将 dist/ 复制到 public/ketcher/
cp -r dist/index.html dist/assets ../public/ketcher/ 2>/dev/null || \
  (mkdir -p ../../public/ketcher && cp -r dist/* ../../public/ketcher/)
cd ../..
```

### 5. 安装 Python 依赖

```bash
pip install rdkit flask flask-cors openpyxl networkx numpy biopython
```

### 6. 数据库初始化

```bash
# Prisma 推送 Schema
bun run db:push
bun run db:generate

# 确保 db/custom.db 存在（已在 05 包中）
```

### 7. 启动服务

```bash
# Next.js 开发服务器 (端口 3000)
bun run dev

# 预测服务 (端口 3003)
cd mini-services/prediction-service
python3 index.py
cd ../..
```

### 8. 访问

```
http://localhost:3000
```

---

## 注意事项

1. **Python 版本**: 需要 3.10+ (推荐 3.12)
2. **Node.js 版本**: 需要 18+ (当前 v24.14.1)
3. **RDKit 安装**: 部分系统可能需要 `conda install -c conda-forge rdkit`
4. **Ketcher**: 已排除 node_modules (约300MB)，需重新 `npm install`
5. **数据库**: `db/custom.db` 已包含全部 51,647 条规则，无需重新导入
6. **Caddyfile**: 需要本地安装 Caddy 才能使用端口转发功能
7. **public/ketcher/**: 已包含预构建的 Ketcher 静态文件，可直接使用
   - 如需重新构建，执行步骤4
