#!/usr/bin/env python3
"""
EZmechanism - 预测服务一键启动脚本

使用方法:
  python start.py          # 启动预测服务
  python start.py --check  # 仅检查依赖，不启动服务
"""

import os
import sys
import subprocess


def check_dependencies():
    """检查所有依赖是否已安装"""
    print("=" * 50)
    print("  EZmechanism 预测服务 - 依赖检查")
    print("=" * 50)
    print()

    missing = []
    version_info = {}

    # 检查 Python 版本
    py_ver = sys.version_info
    py_ok = py_ver >= (3, 9)
    version_info["Python"] = f"{py_ver.major}.{py_ver.minor}.{py_ver.micro}"
    print(f"  {'[OK]' if py_ok else '[FAIL]'} Python {version_info['Python']} {'(需要 >= 3.9)' if not py_ok else ''}")
    if not py_ok:
        missing.append("Python >= 3.9")

    # 检查 pip 包
    packages = {
        "rdkit": "rdkit",
        "flask": "Flask",
        "flask_cors": "flask-cors",
        "openpyxl": "openpyxl",
        "networkx": "networkx",
    }

    for module, display_name in packages.items():
        try:
            mod = __import__(module)
            ver = getattr(mod, "__version__", getattr(mod, "version", ""))
            version_info[display_name] = str(ver) if ver else "OK"
            print(f"  [OK] {display_name} {version_info[display_name]}")
        except ImportError:
            version_info[display_name] = "未安装"
            print(f"  [FAIL] {display_name} 未安装")
            missing.append(display_name)

    # 检查规则文件
    print()
    print("--- 规则文件 ---")
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    rules_files = ["rules_nat_met.xlsx", "rules_test.xlsx"]
    for rf in rules_files:
        path = os.path.join(data_dir, rf)
        if os.path.exists(path):
            size_mb = os.path.getsize(path) / (1024 * 1024)
            print(f"  [OK] {rf} ({size_mb:.1f} MB)")
        else:
            print(f"  [FAIL] {rf} 不存在")
            print(f"        期望路径: {path}")
            missing.append(rf)

    print()
    if missing:
        print("=" * 50)
        print(f"  缺少 {len(missing)} 项依赖，请先安装：")
        print("=" * 50)
        print()

        pip_missing = [m for m in missing if m not in rules_files]
        if pip_missing:
            print("  pip install " + " ".join(pip_missing))

        file_missing = [m for m in missing if m in rules_files]
        if file_missing:
            print(f"  规则文件缺失，请确保 data/ 目录下有: {', '.join(file_missing)}")

        print()
        return False
    else:
        print("=" * 50)
        print("  所有依赖检查通过!")
        print("=" * 50)
        print()
        return True


def start_service():
    """启动 Flask 预测服务"""
    # 切换工作目录到脚本所在目录
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    sys.path.insert(0, script_dir)

    # 启动 index.py
    print("=" * 50)
    print("  启动 M-CSA 预测服务...")
    print("=" * 50)
    print(f"  工作目录: {script_dir}")
    print(f"  Python:   {sys.executable}")
    print(f"  端口:     3003")
    print(f"  访问:     http://localhost:3003")
    print()
    print("  按 Ctrl+C 停止服务")
    print("=" * 50)
    print()

    # 直接导入并启动
    from index import app, load_rules
    load_rules()
    app.run(host="0.0.0.0", port=3003, debug=False, threaded=True)


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    if "--check" in sys.argv or "-c" in sys.argv:
        ok = check_dependencies()
        sys.exit(0 if ok else 1)
        return

    # 检查依赖
    if not check_dependencies():
        print("  依赖检查未通过，请先安装缺少的依赖后重试。")
        sys.exit(1)
        return

    # 启动服务
    start_service()


if __name__ == "__main__":
    main()
