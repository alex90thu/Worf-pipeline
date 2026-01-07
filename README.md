# WORF-Seq Pipeline

WORF-Seq 是一个基于 shell + Python 的全基因组测序（WGS）小型流水线，用于：
- 质控：`fastp`
- 比对：`minimap2`
- BAM 处理：`samtools sort/index`
- 可视化与统计：`scripts/WGSmapping.py`（Matplotlib + NumPy）

该流水线按模块顺序运行：`01_qc.sh → 02_align.sh → 03_process.sh → 04_plot.sh`，由入口脚本 `worf.sh` 调度。

## 目录结构
```
worf.sh
modules/
  01_qc.sh
  02_align.sh
  03_process.sh
  04_plot.sh
references/
scripts/
  WGSmapping.py
```

## 环境准备
建议使用 Conda 创建并激活环境（名称与 `worf.sh` 中一致：`worf-env`）。

```bash
conda env create -f environment.yml
conda activate worf-env
```

如果你只需要 Python 依赖，也可以：

```bash
python -m pip install -r requirements.txt
```

## 外部工具依赖
- fastp（质控）
- minimap2（短读段比对，`-ax sr`）
- samtools（BAM 处理与区域计数）

以上已在 `environment.yml` 中列出；也可以手动安装，例如：

```bash
conda install -c bioconda fastp minimap2 samtools
```

## 输入约定
- 输入目录 `-f <folder>` 必须包含以下原始 FASTQ：
  - `<folder_basename>_raw_1.fq.gz`
  - `<folder_basename>_raw_2.fq.gz`
- `folder_basename` 为输入目录名的最后一段（`basename <folder>`）。

示例：
```
/fdata/sample001/
  sample001_raw_1.fq.gz
  sample001_raw_2.fq.gz
```

## 参考基因组
`modules/02_align.sh` 中的参考路径目前写死为：`/data/guozehua/Worf/references`，并在该目录下优先使用：
- `hg38.mmi`（minimap2 索引）或
- `hg38.fa`（FASTA）

请在该路径准备好参考文件，或根据你的环境修改 `modules/02_align.sh` 的 `REF_DIR`。

创建 minimap2 索引示例：
```bash
minimap2 -d hg38.mmi hg38.fa
```

## 运行方式
入口脚本：`worf.sh`

```bash
bash worf.sh -f <input_folder> -c <chrom> -p <center_bp> \
  [-s step_bp] [-b true|false] [-o outdir] [-r]
```

- `-f`：输入目录路径（包含原始 FASTQ）
- `-c`：染色体名称（如：`chr6`）
- `-p`：中心位置（bp）
- `-s`：步长（默认 `100000` bp，用于全染色体背景分析）
- `-b`：是否进行背景分析（默认 `true`）
- `-o`：输出根目录（默认 `<input_folder>/output`）
- `-r`：强制重新运行（忽略完成标记）

脚本会创建以下目录结构：
```
<outdir>/
  pipeline.log
  01_qc/
  02_align/
  03_bam/
  04_result/
```

## 输出说明
- 质控输出：`01_qc/<sample>_clean_1.fq.gz`, `01_qc/<sample>_clean_2.fq.gz`
- 比对输出：`02_align/<sample>_aligned.sam`
- 处理输出：`03_bam/<sample>_aligned.sorted.bam` 及其 `.bai`
- 结果输出（由 `WGSmapping.py` 生成）：
  - `04_result/<sample>_chromosome_<chr>_step<step>.png`（全染色体覆盖度）
  - `04_result/<sample>_target_region_<chr>_<center>.png`（±50kb精细区域）
  - `04_result/<sample>_worf_seq_summary.txt`（摘要报告，包括中心点附近统计）

## 可视化脚本说明（`scripts/WGSmapping.py`）
- 调用 `samtools view -c` 统计每个 bin 区间的 reads 数
- 使用 Matplotlib 生成覆盖度图；在目标位置处绘制红色虚线并标注坐标
- 需要 Python 包：`matplotlib`, `numpy`

手动运行示例：
```bash
python3 scripts/WGSmapping.py \
  --bam <03_bam/sample_aligned.sorted.bam> \
  --chromosome chr6 \
  --center 32600000 \
  --step 100000 \
  --background true \
  --output <04_result>
```

## 常见问题
- 找不到 `fastp/minimap2/samtools`：请确保 Conda 环境已激活，或手动安装这些工具。
- 找不到原始 FASTQ：确认文件名是否满足约定（`<basename>_raw_1.fq.gz` 等）。
- 参考缺失：在 `modules/02_align.sh` 中修改 `REF_DIR`，指向你的参考文件位置。
- 图像为空或覆盖度为零：确认 BAM 已正确生成并索引，染色体名称与 BAM 头一致（`@SQ SN:<chrom>`）。

## 引用与致谢
- Chen et al., fastp: an ultra-fast all-in-one FASTQ preprocessor.
- Li, Minimap2: pairwise alignment for nucleotide sequences.
- Li et al., The Sequence Alignment/Map format and SAMtools.

---

如需将参考路径做成可配置项或扩展更多分析模块，我可以继续协助完善脚本与文档。