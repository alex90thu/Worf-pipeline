import argparse
from pymongo import MongoClient
from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich.prompt import Prompt, IntPrompt, Confirm
from rich.progress import track
from datetime import datetime
import sys
import re
import os
import csv  # 新增: 用于CSV导出

# 初始化 Rich 控制台
console = Console()

# =======================
# 配置常量
# =======================
EXPORT_DIR = "/data/lulab_commonspace/guozehua/Worf/export"

class WorfDB:
    def __init__(self, host='127.0.0.1', port=30001, db_name='worf_db'):
        self.client = MongoClient(f'mongodb://{host}:{port}/')
        self.db = self.client[db_name]
        self.col = self.db['reads_mapping']

    def get_samples(self):
        """获取所有样品列表"""
        pipeline = [
            {"$group": {
                "_id": "$s",
                "first_doc_id": {"$first": "$_id"}
            }},
            {"$sort": {"_id": 1}}
        ]
        samples = []
        try:
            results = list(self.col.aggregate(pipeline))
            for res in results:
                creation_time = res['first_doc_id'].generation_time.astimezone()
                samples.append({
                    "name": res['_id'],
                    "time": creation_time.strftime("%Y-%m-%d %H:%M:%S")
                })
        except Exception as e:
            console.print(f"[bold red]数据库错误: {e}[/bold red]")
            sys.exit(1)
        return samples

    def get_sample_stats(self, sample_name):
        total_reads = self.col.count_documents({"s": sample_name})
        pipeline = [
            {"$match": {"s": sample_name}},
            {"$group": {
                "_id": "$c",
                "count": {"$sum": 1},
                "avg_read_len": {"$avg": "$rl"},
                "avg_map_len": {"$avg": "$ml"},
                "avg_quality": {"$avg": "$mq"}
            }}
        ]
        chrom_stats = list(self.col.aggregate(pipeline))
        return total_reads, chrom_stats

    def query_reads(self, sample_name, chrom, start=None, end=None):
        query = {"s": sample_name, "c": chrom}
        if end is not None:
            query["st"] = {"$lt": end}
        if start is not None:
            query["ed"] = {"$gt": start}
        
        cursor = self.col.find(query).sort("st", 1)
        return list(cursor)

# =======================
# 辅助函数
# =======================

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]

def is_standard_chrom(chrom_name):
    c = chrom_name.lower()
    if c.startswith('chr'): c = c[3:]
    if c.isdigit(): return True
    if c in ['x', 'y', 'm', 'mt', 'c', 'pt']: return True
    if re.match(r'^[ivx]+$', c): return True
    return False

def parse_region(region_str):
    region_str = region_str.strip().replace(",", "")
    if ":" not in region_str:
        return region_str, None, None
    parts = region_str.split(":", 1)
    chrom = parts[0]
    range_part = parts[1].strip()
    if not range_part: return chrom, None, None
    
    if "-" in range_part:
        r_parts = range_part.split("-")
        start = int(r_parts[0]) if r_parts[0].strip() else None
        end = int(r_parts[1]) if r_parts[1].strip() else None
    else:
        try:
            start = int(range_part)
            end = None
        except:
            return chrom, None, None
    return chrom, start, end

def export_to_csv(reads, sample, chrom, start, end):
    """
    修改点 2 & 3: 导出为 CSV 并保存到固定路径
    """
    # 1. 确保目录存在
    if not os.path.exists(EXPORT_DIR):
        try:
            os.makedirs(EXPORT_DIR, exist_ok=True)
            console.print(f"[dim]已创建导出目录: {EXPORT_DIR}[/dim]")
        except OSError as e:
            console.print(f"[bold red]无法创建目录 {EXPORT_DIR}: {e}[/bold red]")
            return

    # 2. 生成文件名
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    start_str = str(start) if start is not None else "0"
    end_str = str(end) if end is not None else "End"
    
    # 清理文件名中的特殊字符
    safe_chrom = str(chrom).replace("/", "_")
    filename = f"{sample}_{safe_chrom}_{start_str}-{end_str}_{timestamp}.csv"
    full_path = os.path.join(EXPORT_DIR, filename)
    
    # 3. 写入 CSV
    try:
        with open(full_path, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Read Name', 'Chrom', 'Start', 'End', 'Read Len', 'Map Len', 'MapQ', 'CIGAR']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for r in reads:
                writer.writerow({
                    'Read Name': r.get('r'),
                    'Chrom': r.get('c'),
                    'Start': r.get('st'),
                    'End': r.get('ed'),
                    'Read Len': r.get('rl'),
                    'Map Len': r.get('ml'),
                    'MapQ': r.get('mq'),
                    'CIGAR': r.get('cig', '') # 如果没存cig字段这里会为空
                })
                
        console.print(f"[bold green]✅ 导出成功:[/bold green] {full_path}")
    except Exception as e:
        console.print(f"[bold red]写入 CSV 失败: {e}[/bold red]")

# =======================
# 界面逻辑
# =======================

def show_dashboard(db, sample_name):
    # ... (此处代码保持不变，与 v2.0 一致) ...
    # 为了节省篇幅，这里略过 dashboard 显示代码
    # 如果需要完整文件请告诉我
    pass 

def display_reads_table(reads, chrom, start, end, total_count=None):
    title = f"Reads: {chrom}"
    if start or end:
        s_txt = start if start else "0"
        e_txt = end if end else "End"
        title += f":{s_txt}-{e_txt}"
    
    if total_count and total_count > len(reads):
        title += f" (Displaying {len(reads)} of {total_count})"

    rtable = Table(title=title)
    rtable.add_column("Read Name", style="dim")
    rtable.add_column("Start", justify="right")
    rtable.add_column("End", justify="right")
    rtable.add_column("MapLen", justify="right")
    rtable.add_column("MapQ", justify="right")
    
    for r in reads:
        rtable.add_row(
            str(r.get('r')), 
            f"{r.get('st'):,}", 
            f"{r.get('ed'):,}", 
            f"{r.get('ml')}",
            str(r.get('mq'))
        )
    console.print(rtable)

def query_loop(db, sample_name):
    while True:
        console.print("[bold]查询指令:[/bold] [cyan]chr6[/cyan], [cyan]chr6:1000-[/cyan], [cyan]chr6:-5000[/cyan]")
        user_input = Prompt.ask(f"[bold yellow]Query[/bold yellow] ('q' to exit)")
        
        if user_input.lower() == 'q':
            break
            
        chrom, start, end = parse_region(user_input)
        
        if not chrom:
            console.print("[red]无法解析输入！[/red]")
            continue
            
        with console.status(f"[green]Searching {chrom}...[/green]"):
            reads = db.query_reads(sample_name, chrom, start, end)
        
        count = len(reads)
        console.print(f"-> 找到 [bold cyan]{count}[/bold cyan] 条 Reads。")
        
        if count == 0:
            continue

        if count > 20:
            # 先显示前20条
            display_reads_table(reads[:20], chrom, start, end, total_count=count)
            
            action = Prompt.ask(
                "结果过多。选择操作", 
                choices=["show all", "export", "cancel"], 
                default="export"
            )
            
            if action == "export":
                export_to_csv(reads, sample_name, chrom, start, end)
            elif action == "show all":
                display_reads_table(reads, chrom, start, end, total_count=count)
                # 显示全部后，依然给一个保存的机会
                if Confirm.ask("是否保存 CSV?", default=False):
                    export_to_csv(reads, sample_name, chrom, start, end)
        else:
            display_reads_table(reads, chrom, start, end, total_count=count)
            if Confirm.ask("是否保存 CSV?", default=False):
                export_to_csv(reads, sample_name, chrom, start, end)

    console.print("[blue]Bye![/blue]")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=30001)
    args = parser.parse_args()
    
    db = WorfDB(port=args.port)
    
    # 1. Select Sample
    console.rule("[bold blue]Worf Genomics Data Platform[/bold blue]")
    samples = db.get_samples()
    if not samples:
        console.print("[red]无数据[/red]")
        return

    table = Table(title="Available Samples")
    table.add_column("ID", justify="center")
    table.add_column("Name")
    table.add_column("Time")
    for idx, s in enumerate(samples):
        table.add_row(str(idx + 1), s['name'], s['time'])
    console.print(table)
    
    choice = IntPrompt.ask("Select", choices=[str(i+1) for i in range(len(samples))])
    sample_name = samples[choice-1]['name']

    # 2. Dashboard (包含修改后的 dashboard 代码块，为了完整运行需要把之前的 dashboard 代码拷进来)
    # 此处直接重用之前的逻辑，但在 production 代码中确保包含 `show_dashboard` 函数
    # ... (这里需要你保证 show_dashboard 函数在文件里) ...
    
    # 临时插入 show_dashboard 逻辑以保证代码可运行
    console.clear()
    console.rule(f"[bold blue]Sample: {sample_name}[/bold blue]")
    with console.status("[bold yellow]正在计算统计数据...", spinner="earth"):
        total, chrom_stats = db.get_sample_stats(sample_name)
    std_chroms = [c for c in chrom_stats if is_standard_chrom(c['_id'])]
    std_chroms.sort(key=lambda x: natural_sort_key(x['_id']))
    other_count = len(chrom_stats) - len(std_chroms)
    grid = Table.grid(expand=True)
    grid.add_column(justify="center", ratio=1)
    grid.add_column(justify="center", ratio=1)
    grid.add_row(
        Panel(f"[bold cyan]{total:,}[/bold cyan]", title="Total Mapped Reads"),
        Panel(f"[bold cyan]{len(std_chroms)}[/bold cyan] Std + [dim]{other_count} Scaffolds[/dim]", title="Chromosomes")
    )
    console.print(grid)
    ctable = Table(title="Standard Chromosome Statistics")
    ctable.add_column("Chr", style="green")
    ctable.add_column("Reads", justify="right")
    ctable.add_column("Avg Read Len", justify="right")
    ctable.add_column("Avg MapQ", justify="right")
    for c in std_chroms:
        ctable.add_row(c['_id'], f"{c['count']:,}", f"{c['avg_read_len']:.1f}", f"{c['avg_quality']:.1f}")
    console.print(ctable)
    console.print("\n")
    
    # 3. Query Loop
    query_loop(db, sample_name)

if __name__ == "__main__":
    main()