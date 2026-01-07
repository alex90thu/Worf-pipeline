#!/bin/bash
# modules/04_plot.sh

# 【改动】输入从 DIR_BAM 读
BAM_FILE="${DIR_BAM}/${FOLDER_BASENAME}_aligned.sorted.bam"
PY_SCRIPT="${PY_SCRIPT_DIR}/WGSmapping.py"

if [[ ! -f "$BAM_FILE" ]]; then echo "[ERROR] BAM file missing."; exit 1; fi

echo "[MODULE: PLOT] Running Python visualization..."

# 【改动】--output 参数改为 DIR_RESULT
python3 "$PY_SCRIPT" \
    --bam "$BAM_FILE" \
    --chromosome "$CHROMOSOME" \
    --center "$CENTER_POSITION" \
    --step "$STEP_SIZE" \
    --background "$BACKGROUND_ANALYSIS" \
    --output "$DIR_RESULT" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then exit 1; fi