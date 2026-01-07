#!/bin/bash
# worf_pipeline/worf.sh
# 0. 环境准备
conda activate worf-env
# 1. 核心权限设置 (全开放)
umask 000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${SCRIPT_DIR}/modules"
PY_SCRIPT_DIR="${SCRIPT_DIR}/scripts"

# 默认参数
STEP_SIZE=100000
BACKGROUND_ANALYSIS=true
FORCE_RUN=false

usage() {
    echo "Usage: $0 -f <folder> -c <chrom> -p <center> [-s step] [-b background] [-o outdir] [-r force]"
    exit 1
}

while getopts ":f:c:p:s:b:o:r" opt; do
    case $opt in
        f) INPUT_FOLDER="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        p) CENTER_POSITION="$OPTARG" ;;
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

FOLDER_BASENAME=$(basename "$INPUT_FOLDER")

# === 【改动点 1】 修改默认输出目录逻辑 ===
# 如果没有指定 -o，则默认放在 input_folder/output
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="${INPUT_FOLDER}/output"
fi

# === 【改动点 2】 定义并创建分类子目录 ===
# 定义子目录变量并 export 给模块使用
export DIR_QC="${OUTPUT_DIR}/01_qc"
export DIR_ALIGN="${OUTPUT_DIR}/02_align"
export DIR_BAM="${OUTPUT_DIR}/03_bam"
export DIR_RESULT="${OUTPUT_DIR}/04_result"

# 创建目录结构
mkdir -p "$OUTPUT_DIR" "$DIR_QC" "$DIR_ALIGN" "$DIR_BAM" "$DIR_RESULT"

# 日志文件放在根 output 目录下
LOG_FILE="${OUTPUT_DIR}/pipeline.log"

log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1" | tee -a "$LOG_FILE"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" | tee -a "$LOG_FILE"; }

log_info "=== WORF-Seq Pipeline Started ==="
log_info "Input: $INPUT_FOLDER"
log_info "Output Structure Created in: $OUTPUT_DIR"

export INPUT_FOLDER OUTPUT_DIR LOG_FILE SCRIPT_DIR PY_SCRIPT_DIR
export CHROMOSOME CENTER_POSITION STEP_SIZE BACKGROUND_ANALYSIS FOLDER_BASENAME

# 模块调度函数
run_module() {
    local step_name="$1"
    local script_path="${MODULE_DIR}/${step_name}"
    # 标记文件统一放在 output 根目录，方便查看进度
    local checkpoint_file="${OUTPUT_DIR}/.${step_name%.*}.done"

    if [[ "$FORCE_RUN" == "false" && -f "$checkpoint_file" ]]; then
        log_info "检测到完成标记 [${step_name%.*}]，跳过。"
        return 0
    fi

    log_info ">>> Running Step: $step_name"
    if /bin/bash "$script_path"; then
        log_info ">>> Step $step_name Completed."
        touch "$checkpoint_file"
    else
        log_error ">>> Step $step_name FAILED!"
        exit 1
    fi
}

# 按顺序执行
run_module "01_qc.sh"
run_module "02_align.sh"
run_module "03_process.sh"
run_module "04_plot.sh"

log_info "=== Pipeline Successfully Completed ==="
# 最后再次确权，防止意外
chmod -R 777 "$OUTPUT_DIR" 2>/dev/null