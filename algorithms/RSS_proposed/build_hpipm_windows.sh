#!/usr/bin/env bash
# build_hpipm_windows.sh
# 在 Windows 上用 MSYS2 编译 blasfeo 静态库 + hpipm 共享库 (libhpipm.dll)
#
# 前置条件:
#   - 已安装 MSYS2 (默认路径 C:\msys64)
#   - 已通过 MSYS2 UCRT64 环境安装: gcc, make, bc
#       pacboy sync:mman-git ucrt64/toolchain msys/make msys/bc
#
# 用法 (在 PowerShell 中):
#   C:\msys64\usr\bin\env.exe MSYSTEM=UCRT64 /usr/bin/bash -lc "/d/PROJECT/RSS_V2/algorithms/RSS_proposed/build_hpipm_windows.sh"
#
# 或在 MSYS2 UCRT64 终端中:
#   cd /d/PROJECT/RSS_V2
#   bash algorithms/RSS_proposed/build_hpipm_windows.sh

set -e

echo "================================================"
echo "  Build blasfeo + hpipm for Windows (MSYS2)"
echo "================================================"

# 定位脚本所在目录 (即 algorithms/RSS_proposed/), 项目根目录是其上两级
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

BLASFEO_DIR="$PROJECT_ROOT/third_party/blasfeo"
HPIPM_DIR="$PROJECT_ROOT/third_party/hpipm"

echo "PROJECT_ROOT = $PROJECT_ROOT"
echo "BLASFEO_DIR  = $BLASFEO_DIR"
echo "HPIPM_DIR    = $HPIPM_DIR"
echo

# 检查目录存在
if [ ! -d "$BLASFEO_DIR" ]; then
    echo "ERROR: blasfeo 目录不存在: $BLASFEO_DIR"
    echo "       请先执行: git submodule update --init --recursive"
    exit 1
fi
if [ ! -d "$HPIPM_DIR" ]; then
    echo "ERROR: hpipm 目录不存在: $HPIPM_DIR"
    echo "       请先执行: git submodule update --init --recursive"
    exit 1
fi

# 检查工具链
if ! command -v gcc >/dev/null 2>&1; then
    echo "ERROR: gcc 未找到. 请在 MSYS2 UCRT64 环境中运行此脚本."
    echo "       或安装: pacboy sync:mman-git ucrt64/toolchain"
    exit 1
fi
if ! command -v make >/dev/null 2>&1; then
    echo "ERROR: make 未找到. 请安装: pacman -S make"
    exit 1
fi

# 编译 blasfeo 静态库
echo "================================================"
echo "  [1/2] 编译 blasfeo 静态库 (libblasfeo.a)"
echo "================================================"
cd "$BLASFEO_DIR"
make -j4 static_library \
    TARGET=GENERIC \
    USE_C99_MATH=1 \
    EXT_DEP=1 \
    OS=WINDOWS \
    BLASFEO_PATH="$BLASFEO_DIR"

if [ ! -f "$BLASFEO_DIR/lib/libblasfeo.a" ]; then
    echo "ERROR: blasfeo 编译失败, 未生成 libblasfeo.a"
    exit 1
fi
echo "blasfeo 静态库编译完成: $BLASFEO_DIR/lib/libblasfeo.a"
echo

# 编译 hpipm 共享库
echo "================================================"
echo "  [2/2] 编译 hpipm 共享库 (libhpipm.so / libhpipm.dll)"
echo "================================================"
cd "$HPIPM_DIR"
make -j4 shared_library \
    TARGET=GENERIC \
    USE_C99_MATH=1 \
    EXT_DEP=1 \
    OS=WINDOWS \
    BLASFEO_PATH="$BLASFEO_DIR"

if [ ! -f "$HPIPM_DIR/lib/libhpipm.so" ]; then
    echo "ERROR: hpipm 编译失败, 未生成 libhpipm.so"
    exit 1
fi
echo "hpipm 共享库编译完成: $HPIPM_DIR/lib/libhpipm.so"
echo

# Windows: 复制 .so 为 .dll (Python ctypes CDLL('libhpipm.dll') 需要此扩展名)
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OS" == "Windows_NT" ]]; then
    cp -f "$HPIPM_DIR/lib/libhpipm.so" "$HPIPM_DIR/lib/libhpipm.dll"
    echo "已复制 libhpipm.so -> libhpipm.dll (供 Windows Python ctypes 加载)"
fi
echo

echo "================================================"
echo "  编译完成!"
echo "================================================"
echo "产物位置:"
echo "  - $BLASFEO_DIR/lib/libblasfeo.a"
echo "  - $HPIPM_DIR/lib/libhpipm.so"
echo "  - $HPIPM_DIR/lib/libhpipm.dll  (Windows Python 使用)"
echo
echo "Python 加载验证:"
echo "  cd $PROJECT_ROOT"
echo "  python -c \"import sys; sys.path.insert(0,'python'); import hpipm_qp_solver; print('HPIPM_OK=', hpipm_qp_solver._HPIPM_OK)\""
