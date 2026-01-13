#!/bin/bash
# worf_tools.sh
# Description: 独立工具箱，用于手动执行 Worf 的数据库入库(Ingest)与查询(Query)功能
# Author: Guo Zehua

# ================= 配置区域 =================
# 获取当前脚本所在目录
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${BASE_DIR}/scripts"

# MongoDB 配置 (与之前保持一致)
MONGO_HOST="127.0.0.1"
MONGO_PORT=30001

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color
# ===========================================

# 检查 Python 脚本是否存在
if [[ ! -f "${SCRIPTS_DIR}/bam2mongo.py" ]] || [[ ! -f "${SCRIPTS_DIR}/query_worf.py" ]]; then
    echo -e "${RED}[Error] 找不到 Python 脚本！请确保本脚本位于项目根目录，且 scripts/ 文件夹完整。${NC}"
    exit 1
fi

# 检查 MongoDB 容器是否运行 (通过 docker ps)
check_mongo() {
    if ! docker ps | grep -q "worf-mongo"; then
        echo -e "${YELLOW}[Warn] 未检测到名为 'worf-mongo' 的容器运行。${NC}"
        echo -e "尝试连接端口 ${MONGO_PORT}..."
        # 这里不强制退出，因为可能容器名改了，但端口是对的
    fi
}

# 函数: 数据入库
do_ingest() {
    echo -e "${BLUE}=== Worf Data Ingestion Tool ===${NC}"
    
    # 如果没有通过参数传入，则交互式询问
    local bam_file="$1"
    local sample_name="$2"

    if [[ -z "$bam_file" ]]; then
        read -e -p "请输入 BAM 文件路径: " bam_file
    fi

    if [[ ! -f "$bam_file" ]]; then
        echo -e "${RED}[Error] 文件不存在: $bam_file${NC}"
        exit 1
    fi

    if [[ -z "$sample_name" ]]; then
        # 默认使用文件名作为样品名
        local filename=$(basename "$bam_file")
        local default_name="${filename%%.*}" # 去掉后缀
        local today=$(date +%Y%m%d)
        
        read -p "请输入样品名称 (Default: ${today}-${default_name}): " input_name
        sample_name="${input_name:-${today}-${default_name}}"
    fi

    echo -e "${GREEN}正在启动入库程序...${NC}"
    echo -e "  BAM:    $bam_file"
    echo -e "  Sample: $sample_name"
    echo -e "  DB:     $MONGO_HOST:$MONGO_PORT"
    echo "----------------------------------------"

    python "${SCRIPTS_DIR}/bam2mongo.py" \
        --bam "$bam_file" \
        --sample "$sample_name" \
        --host "$MONGO_HOST" \
        --port "$MONGO_PORT"
}

# 函数: 数据查询 (TUI)
do_query() {
    echo -e "${BLUE}=== Worf Interactive Query Interface ===${NC}"
    python "${SCRIPTS_DIR}/query_worf.py" --port "$MONGO_PORT"
}

# 函数: 帮助信息
show_help() {
    echo -e "用法: ./worf_tools.sh [command] [args]"
    echo ""
    echo "Commands:"
    echo "  ingest [bam_path] [sample_id]   将 BAM 文件导入数据库"
    echo "  query                           打开交互式查询面板 (TUI)"
    echo "  help                            显示此帮助"
    echo ""
    echo "Examples:"
    echo "  ./worf_tools.sh query"
    echo "  ./worf_tools.sh ingest ./data/test.bam MySample01"
}

# ================= 主逻辑 =================
check_mongo

case "$1" in
    ingest)
        do_ingest "$2" "$3"
        ;;
    query)
        do_query
        ;;
    *)
        if [[ -z "$1" ]]; then
            # 如果没输参数，给个简单的选择菜单
            echo "请选择操作模式:"
            echo "1) Query  (数据查询)"
            echo "2) Ingest (数据入库)"
            read -p "Select [1-2]: " choice
            case "$choice" in
                1) do_query ;;
                2) do_ingest ;;
                *) echo -e "${RED}无效选择${NC}"; exit 1 ;;
            esac
        else
            show_help
            exit 1
        fi
        ;;
esac