import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from datetime import datetime
import argparse
import re

def get_counts(bam_file, chrom, start, end, bin_size):
    """使用samtools计算指定区间内每个bin的符合条件的reads数"""
    bins = range(start, end, bin_size)
    counts = []
    
    # 检查samtools是否可用
    try:
        subprocess.run(['samtools', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        raise FileNotFoundError("samtools not found in PATH") from e
    
    for b_start in bins:
        b_end = b_start + bin_size
        bin_count = 0
        try:
            # 使用samtools view获取指定区域的reads
            cmd = [
                'samtools', 'view', '-c', 
                bam_file, 
                f'{chrom}:{b_start}-{b_end-1}'
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # 移除空白字符并转换为整数
            bin_count = int(result.stdout.strip())
        except (subprocess.CalledProcessError, ValueError):
            bin_count = 0
            
        counts.append(bin_count)
        
    return list(bins), counts

def calculate_flanking_stats(bam_file, chrom, center, radii=[200, 500, 1000]):
    """
    计算中心点上下游特定范围内的Reads总数
    radii: 距离中心点的半径列表 (e.g., [200] 意味着范围是 center-200 到 center+200)
    """
    stats = {}
    print(f"\n[INFO] 正在统计中心点 ({chrom}:{center}) 附近的 reads 数...")
    
    for r in radii:
        start = max(0, center - r)
        end = center + r
        region_str = f"{chrom}:{start}-{end}"
        try:
            cmd = ['samtools', 'view', '-c', bam_file, region_str]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            count = int(result.stdout.strip())
            stats[r] = count
            print(f"  - +/- {r}bp (Range: {start}-{end}): {count:,} reads")
        except subprocess.CalledProcessError:
            print(f"  [WARN] 无法计算区域 {region_str}")
            stats[r] = "N/A"
            
    return stats

def plot_data(bins, counts, chrom, bin_size, title, filename, target_pos=None):
    """绘图并保存 (改进版配色)"""
    
    plt.figure(figsize=(12, 5))
    
    # 转换为Mb单位
    x_centers = (np.array(bins) + bin_size / 2.0) / 1e6
    counts_arr = np.array(counts)

    # --- [改进点1] 配色方案修改 ---
    # 放弃原本不可见的 RdYlBu_r，改为高对比度的纯色 'SteelBlue'
    # 这种颜色在白色背景下清晰，且符合学术出版标准
    bar_color = 'steelblue' 
    
    # 绘制线条
    # 使用 vlines 绘制，颜色统一
    plt.vlines(x_centers, 0, counts_arr, color=bar_color, linewidth=0.9, alpha=0.9)
    
    # 对于计数为0的区域，如果不画任何东西可能会导致断裂感，
    # 但在覆盖度图中，留白通常就是表示0。
    # 如果为了美观需要底线，可以取消下面这行的注释：
    # plt.axhline(0, color='gray', linewidth=0.5)

    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel(f"Chromosome {chrom} Position (Mb)", fontsize=12)
    plt.ylabel("Read Counts", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.3)

    # 如果给定目标位置，则绘制垂直虚线并标注原始坐标
    if target_pos is not None:
        x_target_mb = target_pos / 1e6
        ymax = counts_arr.max() if counts_arr.size else 1
        
        # 绘制目标位置红线
        plt.axvline(x=x_target_mb, color='#e41a1c', linestyle='--', linewidth=1.5, alpha=0.8)
        
        # 标注文字
        plt.text(x_target_mb, ymax * 0.95, f" {int(target_pos):,}", rotation=90,
                 va='top', ha='right', color='#e41a1c', fontsize=10, fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#e41a1c", alpha=0.8))

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"✅ 成功生成图像: {filename}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='WGSmapping: WGS background and target enrichment plotting')
    parser.add_argument('--bam', required=True, help='Input BAM file path')
    parser.add_argument('--chromosome', required=True, help='Target chromosome (e.g., chr6)')
    parser.add_argument('--center', type=int, required=True, help='Center position (bp)')
    parser.add_argument('--step', type=int, default=100000, help='Step size for genome-wide analysis (bp)')
    parser.add_argument('--background', type=str, default='true', help='Perform background analysis (true/false)')
    parser.add_argument('--output', required=True, help='Output directory for plots')

    args = parser.parse_args()

    print("[INFO] === WORF-Seq 染色体比对分析 (v2.0) ===")
    
    # 1. 检查BAM文件
    bam_path = args.bam
    if not os.path.exists(bam_path):
        print(f"[ERROR] BAM文件不存在: {bam_path}")
        return
    
    # 2. 输出目录设置
    out_dir = args.output
    os.makedirs(out_dir, exist_ok=True)
    bam_basename = os.path.splitext(os.path.basename(bam_path))[0]
    sample_prefix = re.sub(r'(_aligned_minimap)?(\.sorted|_sorted)?$', '', bam_basename)

    target_chrom = args.chromosome
    target_pos = args.center
    wgs_bin = args.step
    do_background = args.background.lower() in ['true', 'yes', '1', 'on']

    # 获取染色体长度
    try:
        cmd = ['samtools', 'view', '-H', bam_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        chrom_length = None
        for line in result.stdout.split('\n'):
            if line.startswith('@SQ') and f'SN:{target_chrom}' in line:
                for part in line.split('\t'):
                    if part.startswith('LN:'):
                        chrom_length = int(part[3:])
                        break
                break
        
        if chrom_length is None:
            print(f"[ERROR] 无法获取染色体 {target_chrom} 的长度")
            return
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] 获取染色体长度失败: {e}")
        return

    generated_files = []

    # --- [改进点2] 计算中心区域统计信息 ---
    # 定义需要统计的半径范围
    stat_radii = [200, 500, 1000]
    flanking_stats = calculate_flanking_stats(bam_path, target_chrom, target_pos, stat_radii)

    # --- 执行全长分析 ---
    if do_background:
        print(f"\n[INFO] [1/2] 正在分析 {target_chrom} 全长背景...")
        try:
            wgs_bins, wgs_counts = get_counts(bam_path, target_chrom, 0, chrom_length, wgs_bin)
            wgs_fname = os.path.join(out_dir, f"{sample_prefix}_chromosome_{target_chrom}_step{wgs_bin}.png")
            plot_data(wgs_bins, wgs_counts, target_chrom, wgs_bin,
                     f"WORF-Seq Chromosome Coverage\\n{target_chrom} (Step: {wgs_bin:,} bp)", 
                     wgs_fname, target_pos=target_pos)
            generated_files.append(wgs_fname)
        except Exception as e:
            print(f"[ERROR] 全染色体分析失败: {e}")
    else:
        print("[INFO] 跳过全染色体分析")

    # --- 执行精细分析 ---
    micro_bin = 500
    micro_start = max(0, target_pos - 50000)
    micro_end = min(chrom_length, target_pos + 50000)

    print(f"[INFO] [2/2] 正在分析目标区域 (+/- 50kb)...")
    try:
        m_bins, m_counts = get_counts(bam_path, target_chrom, micro_start, micro_end, micro_bin)
        target_fname = os.path.join(out_dir, f"{sample_prefix}_target_region_{target_chrom}_{target_pos}.png")
        plot_data(m_bins, m_counts, target_chrom, micro_bin,
                 f"WORF-Seq Target Region Coverage\\n{target_chrom}:{micro_start:,}-{micro_end:,}", 
                 target_fname, target_pos=target_pos)
        generated_files.append(target_fname)
    except Exception as e:
        print(f"[ERROR] 目标区域分析失败: {e}")
    
    # --- 生成摘要报告 (包含统计信息) ---
    summary_fname = os.path.join(out_dir, f"{sample_prefix}_worf_seq_summary.txt")
    try:
        with open(summary_fname, 'w') as f:
            f.write("WORF-Seq Analysis Summary Report\\n")
            f.write("=" * 40 + "\\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
            f.write(f"BAM File: {bam_path}\\n")
            f.write(f"Target: {target_chrom}:{target_pos:,}\\n\\n")
            
            f.write("Target Region Statistics (Cumulative Reads):\\n")
            f.write("-" * 40 + "\\n")
            for r in stat_radii:
                count = flanking_stats.get(r, "N/A")
                range_str = f"{target_chrom}:{max(0, target_pos-r)}-{target_pos+r}"
                f.write(f"  +/- {r:<4} bp (Total span {2*r:<4} bp): {count:>8} reads  [{range_str}]\\n")
            f.write("-" * 40 + "\\n\\n")

            f.write("Generated Files:\\n")
            for file in generated_files:
                f.write(f"- {file}\\n")
        generated_files.append(summary_fname)
        print(f"[INFO] 摘要报告已保存: {summary_fname}")
    except Exception as e:
        print(f"[ERROR] 生成摘要报告失败: {e}")

    print(f"\n[SUCCESS] 分析完成，生成 {len(generated_files)} 个文件。")

if __name__ == "__main__":
    main()