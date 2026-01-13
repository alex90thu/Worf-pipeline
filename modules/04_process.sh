#!/bin/bash
# modules/03_process.sh

# 【改动】输入从 DIR_ALIGN 读，输出到 DIR_BAM
SAM_FILE="${DIR_ALIGN}/${FOLDER_BASENAME}_aligned.sam"
BAM_FILE="${DIR_BAM}/${FOLDER_BASENAME}_aligned.sorted.bam"

if [[ ! -f "$SAM_FILE" ]]; then echo "[ERROR] SAM file missing."; exit 1; fi

echo "[MODULE: PROCESS] Sorting BAM..."
samtools sort -@ 8 -o "$BAM_FILE" "$SAM_FILE" 2>&1 | tee -a "$LOG_FILE"

if [[ ${PIPESTATUS[0]} -eq 0 ]]; then
    echo "[MODULE: PROCESS] Indexing BAM..."
    samtools index "$BAM_FILE" 2>&1 | tee -a "$LOG_FILE"
else
    echo "[ERROR] BAM Sort failed."; exit 1
fi