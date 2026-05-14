#!/usr/bin/env python3
"""
EZmechanism - 一键启动脚本（前端 + 预测服务）

使用方法:
  python start_all.py          # 同时启动前端和预测服务
  python start_all.py --frontend   # 只启动前端
  python start_all.py --backend    # 只启动预测服务
"""

import os
import sys
import subprocess
import time
import platform
import signal


# 获取项目根目录（src 的上一级）
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
PREDICTION_DIR = SCRIPT_DIR

processes = []


def is_windows():
    return platform.system() == "Windows"


def cleanup(signum=None, frame=None):
    """停止所有子进程"""
    print("\n\n  正在停止所有服务...")
    for p in processes:
        try:
            if is_windows():
                p.terminate()
            else:
                os.killpg(os.getpgid(p.pid), signal.SIGTERM)
        except Exception:
            pass
    print("  所有服务已停止。")
    sys.exit(0)


signal.signal(signal.SIGINT, cleanup)
signal.signal(signal.SIGTERM, cleanup)


def check_node():
    """检查 Node.js 是否安装"""
    try:
        result = subprocess.run(
            ["node", "--version"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            ver = result.stdout.strip()
            print(f"  [OK] Node.js {ver}")
            return True
    except Exception:
        pass
    print("  [FAIL] Node.js 未安装，前端无法启动")
    print("         请先安装: https://nodejs.org")
    return False


def check_npm():
    """检查 npm/bun 是否可用"""
    for cmd in ["bun", "npm"]:
        try:
            result = subprocess.run(
                [cmd, "--version"],
                capture_output=True, text=True, timeout=5
            )
            if result.returncode == 0:
                ver = result.stdout.strip()
                print(f"  [OK] {cmd} {ver}")
                return cmd
        except Exception:
            pass
    return None


def check_database():
    """检查数据库是否已初始化"""
    db_path = os.path.join(PROJECT_ROOT, "db", "mcsa_rules.db")
    if os.path.exists(db_path):
        print(f"  [OK] 数据库已存在")
        return True
    else:
        print(f"  [WARN] 数据库不存在: {db_path}")
        print(f"         前端首次启动时会自动创建")
        return True


def check_python_deps():
    """检查 Python 依赖"""
    missing = []
    packages = {"rdkit": "rdkit", "flask": "Flask", "flask_cors": "flask-cors",
                "openpyxl": "openpyxl", "networkx": "networkx"}
    for module, name in packages.items():
        try:
            __import__(module)
            print(f"  [OK] {name}")
        except ImportError:
            print(f"  [FAIL] {name}")
            missing.append(name)

    if missing:
        print(f"\n  缺少 Python 依赖，正在安装: {', '.join(missing)}")
        try:
            subprocess.run(
                [sys.executable, "-m", "pip", "install"] + missing,
                check=True
            )
            print("  安装完成!")
            return True
        except Exception as e:
            print(f"  安装失败: {e}")
            return False
    return True


def check_rules_files():
    """检查规则文件"""
    data_dir = os.path.join(PREDICTION_DIR, "data")
    for f in ["rules_nat_met.xlsx"]:
        path = os.path.join(data_dir, f)
        if os.path.exists(path):
            size_mb = os.path.getsize(path) / (1024 * 1024)
            print(f"  [OK] {f} ({size_mb:.1f} MB)")
        else:
            print(f"  [WARN] {f} 不存在，预测功能将不可用")


def start_frontend():
    """启动 Next.js 前端"""
    print("\n" + "=" * 50)
    print("  启动前端服务 (Next.js)...")
    print("=" * 50)

    npm_cmd = check_npm()
    if not npm_cmd:
        print("  [FAIL] 未找到 npm 或 bun，跳过前端")
        return

    os.chdir(PROJECT_ROOT)

    # 检查 node_modules
    if not os.path.exists(os.path.join(PROJECT_ROOT, "node_modules")):
        print("  node_modules 不存在，正在安装依赖...")
        subprocess.run([npm_cmd, "install"], cwd=PROJECT_ROOT)

    print(f"  启动命令: {npm_cmd} run dev")
    print(f"  工作目录: {PROJECT_ROOT}")
    print(f"  访问地址: http://localhost:3000")

    if is_windows():
        p = subprocess.Popen(
            [npm_cmd, "run", "dev"],
            cwd=PROJECT_ROOT,
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP,
        )
    else:
        p = subprocess.Popen(
            [npm_cmd, "run", "dev"],
            cwd=PROJECT_ROOT,
            preexec_fn=os.setsid,
        )
    processes.append(p)
    print(f"  前端已启动 (PID: {p.pid})")


def start_backend():
    """启动预测服务"""
    print("\n" + "=" * 50)
    print("  启动预测服务 (Flask)...")
    print("=" * 50)

    os.chdir(PREDICTION_DIR)
    sys.path.insert(0, PREDICTION_DIR)

    print(f"  启动命令: python index.py")
    print(f"  工作目录: {PREDICTION_DIR}")
    print(f"  服务地址: http://localhost:3003")

    if is_windows():
        p = subprocess.Popen(
            [sys.executable, "index.py"],
            cwd=PREDICTION_DIR,
            creationflags=subprocess.CREATE_NEW_PROCESS_GROUP,
        )
    else:
        p = subprocess.Popen(
            [sys.executable, "index.py"],
            cwd=PREDICTION_DIR,
            preexec_fn=os.setsid,
        )
    processes.append(p)
    print(f"  预测服务已启动 (PID: {p.pid})")


def main():
    only_frontend = "--frontend" in sys.argv
    only_backend = "--backend" in sys.argv

    print()
    print("=" * 50)
    print("  EZmechanism 一键启动")
    print("=" * 50)
    print(f"  项目目录: {PROJECT_ROOT}")
    print(f"  预测服务: {PREDICTION_DIR}")
    print()

    # ---- 检查依赖 ----
    if not only_backend:
        print("--- 前端检查 ---")
        check_node()
        check_npm()
        check_database()

    if not only_frontend:
        print("\n--- 后端检查 ---")
        check_python_deps()
        check_rules_files()

    # ---- 启动服务 ----
    if only_frontend:
        start_frontend()
    elif only_backend:
        start_backend()
    else:
        start_backend()
        time.sleep(1)
        start_frontend()

    # ---- 显示摘要 ----
    print("\n" + "=" * 50)
    print("  所有服务启动完毕!")
    print("=" * 50)
    if not only_backend:
        print("  前端页面:  http://localhost:3000")
    if not only_frontend:
        print("  预测服务:  http://localhost:3003/api/health")
        print("  规则数量:  51647+")
    print()
    print("  按 Ctrl+C 停止所有服务")
    print("=" * 50)
    print()

    # 等待所有进程
    while True:
        for p in processes[:]:
            ret = p.poll()
            if ret is not None:
                print(f"  [WARN] 进程 PID {p.pid} 已退出 (code: {ret})")
                processes.remove(p)
        if not processes:
            print("  所有服务已停止。")
            break
        time.sleep(2)


if __name__ == "__main__":
    main()
