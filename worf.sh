#!/bin/bash
# worf_pipeline/worf.sh
umask 000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULE_DIR="${SCRIPT_DIR}/modules"
PY_SCRIPT_DIR="${SCRIPT_DIR}/scripts"

# === [修改点] 默认参数增加 WINDOW_SIZE ===
STEP_SIZE=100000
WINDOW_SIZE=10000  # 默认精细窗口半径 (10kb)
BACKGROUND_ANALYSIS=true
FORCE_RUN=false

usage() {
    # === [修改点] 帮助信息 ===
    echo "Usage: $0 -f <folder> -c <chrom> -p <center> [-w window] [-s step] [-b background] [-o outdir]"
    echo "  -w  Target window radius (bp) for linear pile-up plot (default: 10000)"
    exit 1
}

# === [修改点] 获取 -w 参数 ===
while getopts ":f:c:p:w:s:b:o:r" opt; do
    case $opt in
        f) INPUT_FOLDER="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        p) CENTER_POSITION="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;  # 获取窗口大小
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

# ... (中间部分保持不变: 目录创建、日志设置等) ...
# 注意：确保导出 WINDOW_SIZE 变量

FOLDER_BASENAME=$(basename "$INPUT_FOLDER")
if [[ -z "$OUTPUT_DIR" ]]; then OUTPUT_DIR="${INPUT_FOLDER}/output"; fi

export DIR_QC="${OUTPUT_DIR}/01_qc"
export DIR_ALIGN="${OUTPUT_DIR}/02_align"
export DIR_BAM="${OUTPUT_DIR}/03_bam"
export DIR_RESULT="${OUTPUT_DIR}/04_result"
mkdir -p "$OUTPUT_DIR" "$DIR_QC" "$DIR_ALIGN" "$DIR_BAM" "$DIR_RESULT"

LOG_FILE="${OUTPUT_DIR}/pipeline.log"
log_info() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1" | tee -a "$LOG_FILE"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" | tee -a "$LOG_FILE"; }

# === [修改点] 导出 WINDOW_SIZE ===
export INPUT_FOLDER OUTPUT_DIR LOG_FILE SCRIPT_DIR PY_SCRIPT_DIR
export CHROMOSOME CENTER_POSITION STEP_SIZE WINDOW_SIZE BACKGROUND_ANALYSIS FOLDER_BASENAME

run_module() {
    local step_name="$1"
    local script_path="${MODULE_DIR}/${step_name}"
    local checkpoint_file="${OUTPUT_DIR}/.${step_name%.*}.done"

    # 对 04_plot.sh 不进行 checkpoint 跳过与创建，保证可重复绘图
    if [[ "$step_name" != "04_plot.sh" && "$FORCE_RUN" == "false" && -f "$checkpoint_file" ]]; then
        log_info "Skipping $step_name (Checkpointed)"
        return 0
    fi
    log_info "Running $step_name..."
    if /bin/bash "$script_path"; then
        log_info "$step_name Completed."
        if [[ "$step_name" != "04_plot.sh" ]]; then
            touch "$checkpoint_file"
        fi
    else
        log_error "$step_name Failed."
        exit 1
    fi
}

# 执行
run_module "01_qc.sh"
run_module "02_align.sh"
run_module "03_process.sh"
run_module "04_plot.sh"

chmod -R 777 "$OUTPUT_DIR" 2>/dev/null
log_info "Pipeline Done."