#!/bin/bash
# worf_pipeline/worf.sh
umask 000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${SCRIPT_DIR}/modules"
PY_SCRIPT_DIR="${SCRIPT_DIR}/scripts"

# === 默认参数 ===
STEP_SIZE=100000
WINDOW_SIZE=10000  # 默认精细窗口半径 (10kb)
BACKGROUND_ANALYSIS=true
FORCE_RUN=false

usage() {
    echo "Usage: $0 -f <folder> -c <chrom> -p <center> [-w window] [-s step] [-b background] [-o outdir]"
    echo "  -w  Target window radius (bp) for linear pile-up plot (default: 10000)"
    echo "  -r  Force re-run (overwrite existing checkpoints)"
    exit 1
}

# === 参数解析 ===
while getopts ":f:c:p:w:s:b:o:r" opt; do
    case $opt in
        f) INPUT_FOLDER="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        p) CENTER_POSITION="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;
        s) STEP_SIZE="$OPTARG" ;;
        b) BACKGROUND_ANALYSIS="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        r) FORCE_RUN=true ;;
        \?) usage ;;
        :) usage ;;
    esac
done

if [[ -z "$INPUT_FOLDER" || -z "$CHROMOSOME" || -z "$CENTER_POSITION" ]]; then
    echo "[ERROR] Missing required arguments."
    usage
fi

# === 核心变量定义 ===
FOLDER_BASENAME=$(basename "$INPUT_FOLDER")
# [NEW] 定义统一的 Sample ID，用于数据库索引
SAMPLE_ID="${FOLDER_BASENAME}"

if [[ -z "$OUTPUT_DIR" ]]; then OUTPUT_DIR="${INPUT_FOLDER}/output"; fi

# 定义子目录结构
export DIR_QC="${OUTPUT_DIR}/01_qc"
export DIR_ALIGN="${OUTPUT_DIR}/02_align"
# export DIR_BAM="${OUTPUT_DIR}/03_bam" # (旧结构，暂保留以防兼容)
export DIR_RESULT="${OUTPUT_DIR}/04_result"

# 创建目录
mkdir -p "$OUTPUT_DIR" "$DIR_QC" "$DIR_ALIGN" "$DIR_RESULT"

# 日志设置
LOG_FILE="${OUTPUT_DIR}/pipeline.log"
log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1" | tee -a "$LOG_FILE"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" | tee -a "$LOG_FILE"; }

# === [重要] 导出所有环境变量供子模块使用 ===
export INPUT_FOLDER OUTPUT_DIR LOG_FILE SCRIPT_DIR PY_SCRIPT_DIR
export CHROMOSOME CENTER_POSITION STEP_SIZE WINDOW_SIZE BACKGROUND_ANALYSIS 
export FOLDER_BASENAME SAMPLE_ID

# === 模块执行函数 ===
run_module() {
    local step_name="$1"
    local script_path="${MODULE_DIR}/${step_name}"
    local checkpoint_file="${OUTPUT_DIR}/.${step_name%.*}.done"

    # 检查脚本是否存在
    if [[ ! -f "$script_path" ]]; then
        log_error "Module script not found: $script_path"
        exit 1
    fi

    # Checkpoint 检查 (04_plot 总是运行，不跳过)
    if [[ "$step_name" != "04_plot.sh" && "$FORCE_RUN" == "false" && -f "$checkpoint_file" ]]; then
        log_info "Skipping $step_name (Checkpointed)"
        return 0
    fi

    log_info "Running $step_name for Sample: $SAMPLE_ID..."
    
    # [修改点] 执行时传入 SAMPLE_ID 作为 $1 参数
    # 这样 modules/03_db_ingest.sh 里的 SAMPLE_ID=$1 就能正确获取
    if /bin/bash "$script_path" "$SAMPLE_ID"; then
        log_info "$step_name Completed."
        # 成功后创建 checkpoint
        if [[ "$step_name" != "04_plot.sh" ]]; then
            touch "$checkpoint_file"
        fi
    else
        log_error "$step_name Failed."
        exit 1
    fi
}

# === Pipeline 执行流 ===

# Step 1: 质控
run_module "01_qc.sh"

# Step 2: 比对 (生成 BAM)
run_module "02_align.sh"

# Step 3: [NEW] 数据入库 (BAM -> MongoDB)
# 替代了旧的 03_process.sh
run_module "03_db_ingest.sh"

# Step 4: 绘图 (针对参数指定的区域进行可视化)
run_module "04_plot.sh"

# 权限修正 (防止 Docker 生成的文件锁死)
chmod -R 777 "$OUTPUT_DIR" 2>/dev/null
log_info "Pipeline Done."