#!/bin/bash
# modules/01_qc.sh

RAW_R1="${INPUT_FOLDER}/${FOLDER_BASENAME}_raw_1.fq.gz"
RAW_R2="${INPUT_FOLDER}/${FOLDER_BASENAME}_raw_2.fq.gz"

# 【改动】输出路径改为 DIR_QC
CLEAN_R1="${DIR_QC}/${FOLDER_BASENAME}_clean_1.fq.gz"
CLEAN_R2="${DIR_QC}/${FOLDER_BASENAME}_clean_2.fq.gz"

echo "[MODULE: QC] Input: $RAW_R1"
echo "[MODULE: QC] Output Dir: $DIR_QC"

if [[ ! -f "$RAW_R1" || ! -f "$RAW_R2" ]]; then
    echo "[ERROR] Raw files not found!"
    exit 1
fi

if command -v fastp >/dev/null 2>&1; then
    fastp -i "$RAW_R1" -I "$RAW_R2" -o "$CLEAN_R1" -O "$CLEAN_R2" -Q -L 2>&1 | tee -a "$LOG_FILE"
    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then exit 1; fi
else
    echo "[WARN] fastp missing, using raw files." | tee -a "$LOG_FILE"
    ln -sf "$RAW_R1" "$CLEAN_R1"
    ln -sf "$RAW_R2" "$CLEAN_R2"
fi