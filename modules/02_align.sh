#!/bin/bash
# modules/02_align.sh

# 允许通过环境变量覆盖参考目录，未设置则使用默认路径
REF_DIR="${REF_DIR:-/data/lulab_commonspace/guozehua/Worf/references}"
HG38_MMI="${REF_DIR}/hg38.mmi"
HG38_FA="${REF_DIR}/hg38.fa"

# 【改动】输入从 DIR_QC 读，输出到 DIR_ALIGN
CLEAN_R1="${DIR_QC}/${FOLDER_BASENAME}_clean_1.fq.gz"
CLEAN_R2="${DIR_QC}/${FOLDER_BASENAME}_clean_2.fq.gz"
SAM_FILE="${DIR_ALIGN}/${FOLDER_BASENAME}_aligned.sam"

if [[ -f "$HG38_MMI" ]]; then REF_FILE="$HG38_MMI";
elif [[ -f "$HG38_FA" ]]; then REF_FILE="$HG38_FA";
else echo "[ERROR] Reference missing in $REF_DIR"; exit 1; fi

if [[ ! -f "$CLEAN_R1" ]]; then echo "[ERROR] Clean data missing."; exit 1; fi

echo "[MODULE: ALIGN] Minimap2 running..."
echo "[MODULE: ALIGN] Using reference dir: $REF_DIR"
minimap2 -ax sr -t 8 "$REF_FILE" "$CLEAN_R1" "$CLEAN_R2" > "$SAM_FILE" 2>> "$LOG_FILE"

if [[ $? -ne 0 ]]; then echo "[ERROR] Alignment failed."; exit 1; fi