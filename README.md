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
建议使用 Conda 创建并激活环境（名称与 `worf.sh` 一致：`worf_env`）。

```bash
conda env create -f environment.yml
conda activate worf_env
```

如果你只需要 Python 依赖，也可以：

```bash
python -m pip install -r requirements.txt
```

## 外部工具依赖
- fastp（质控）
- minimap2（短读段比对，`-ax sr`）
- samtools（BAM 处理与区域计数）
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
建议使用 Conda 创建并激活环境（默认环境名：`worf_env`）。

```bash
conda env create -f environment.yml
conda activate worf_env
```

如果只需要 Python 包，也可以使用 pip：

```bash
python -m pip install -r requirements.txt
```

## 外部工具依赖
- `fastp`（质控）
- `minimap2`（比对）
- `samtools`（BAM 处理与区域计数）

这些工具已在 [environment.yml](environment.yml) 中列出。若使用 Conda 手动安装：

```bash
conda install -c bioconda fastp minimap2 samtools
```

注意：运行时需要确保 `samtools` 在 `PATH` 中或已激活包含 samtools 的环境。

## 运行方式
入口脚本：`worf.sh`

```bash
bash worf.sh -f <input_folder> -c <chrom> -p <center_bp> \
  [-w window_bp] [-s step_bp] [-b true|false] [-o outdir] [-r]
```

- `-f`：输入目录（包含原始 FASTQ）
- `-c`：染色体（例如 `chr6` 或 `6`，需与 BAM 头的 `@SQ SN:` 一致）
- `-p`：中心位置（bp）
- `-w`：目标窗口半径（bp，默认 `10000`）
- `-s`：步长（bp，默认 `100000`，用于全染色体背景统计）
- `-b`：是否进行背景分析（默认 `true`）
- `-o`：输出目录（默认 `<input_folder>/output`）
- `-r`：强制重跑（忽略检查点）

示例：

```bash
bash worf.sh -f /data/.../UDI001 -c chr6 -p 31429000 -w 10000
```

### Checkpoints 与重复绘图
- 前三步（`01_qc.sh`, `02_align.sh`, `03_process.sh`）仍使用 checkpoint 机制（`.01_qc.done` 等）。
- `04_plot.sh` 不再创建/检查 checkpoint，允许使用不同 `-c/-p/-w` 参数反复绘图。

### 绘图输出命名
- 由 `scripts/WGSmapping.py` 生成的图像文件名包含 `-c`（染色体）、`-p`（中心位置）、`-w`（窗口半径）与时间戳，便于区分多次绘图。
- 命名示例：
  - `SAMPLE_cchr6_p31249000_w10000_20260108_142233_hist.png`
  - `SAMPLE_cchr6_p31249000_w10000_20260108_142233_pileup.png`

脚本会创建以下目录（在 `<outdir>` 下）：

```
pipeline.log
01_qc/
02_align/
03_bam/
04_result/
```

## 可视化脚本说明（`scripts/WGSmapping.py`）
- 使用 `samtools` 统计区域内 Reads 并生成覆盖度图与堆叠图。
- 脚本优先使用 `BrokenBarHCollection` 提升绘图性能；若当前 matplotlib 版本缺失该类，脚本会自动回退到 `PatchCollection` + `Rectangle` 的实现以保证兼容性。

手动运行示例：

```bash
python3 scripts/WGSmapping.py \
  --bam 03_bam/sample_aligned.sorted.bam \
  --chromosome chr6 \
  --center 31429000 \
  --step 100000 \
  --background true \
  --output 04_result
```

## 常见问题
- 找不到 `fastp/minimap2/samtools`：请确保 Conda 环境已激活，或手动安装这些工具。
- 找不到原始 FASTQ：确认文件名是否满足约定（`<basename>_raw_1.fq.gz` 等）。
- 参考缺失：在 `modules/02_align.sh` 中修改 `REF_DIR`，指向你的参考文件位置。
- 图像为空或覆盖度为零：确认 BAM 已正确生成并索引，染色体名称与 BAM 头一致（`@SQ SN:<chrom>`）。

## 致谢
- fastp, minimap2, samtools; Matplotlib 与 NumPy 用于绘图与数值计算。

如需我把参考路径参数化或添加更多分析模块，我可以继续完善脚本与文档。