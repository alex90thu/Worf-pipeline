"""fetch_genes.py

小工具：通过 NCBI E-utilities 查询基因坐标并以命令行工具形式暴露。

用法示例：
  python fetch_genes.py -g TP53
  echo -e "TP53\nBRCA1" | python fetch_genes.py --format json
  python fetch_genes.py -f genes.txt --format tsv
"""

from __future__ import annotations

import argparse
import sys
import json
import csv
import requests
from typing import Optional, Tuple, Dict, Any, List


def log_gene_lookup(msg: str) -> None:
    """简单日志函数，写到 stderr。"""
    sys.stderr.write(f"[fetch_genes] {msg}\n")


def fetch_gene_coordinates(gene_symbol: str, organism: str = "Homo sapiens") -> Tuple[Optional[Dict[str, Any]], str]:
    """通过NCBI E-utilities查询基因坐标，返回字典或 None 和信息字符串。"""
    if not gene_symbol:
        return None, "请输入基因名称"

    try:
        query_term = f"{gene_symbol}[gene] AND {organism}[organism]"
        log_gene_lookup(f"esearch term='{query_term}'")

        esearch_params = {
            "db": "gene",
            "term": query_term,
            "retmode": "json",
            "retmax": 5,
        }
        esearch_resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
            params=esearch_params,
            timeout=10,
        )
        log_gene_lookup(f"esearch status={esearch_resp.status_code} url={esearch_resp.url}")
        esearch_resp.raise_for_status()
        esearch_data = esearch_resp.json()
        id_list = esearch_data.get("esearchresult", {}).get("idlist", [])
        log_gene_lookup(f"esearch idlist={id_list}")
        if not id_list:
            return None, "未找到匹配的基因，请输入官方基因符号（如 TP53, HLA-C）"

        gene_id = id_list[0]
        esummary_params = {"db": "gene", "id": gene_id, "retmode": "json"}
        esummary_resp = requests.get(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
            params=esummary_params,
            timeout=10,
        )
        log_gene_lookup(f"esummary status={esummary_resp.status_code} url={esummary_resp.url}")
        esummary_resp.raise_for_status()
        esummary_data = esummary_resp.json()

        docsum = esummary_data.get("result", {}).get(str(gene_id), {})
        if not docsum and "DocumentSummarySet" in esummary_data:
            summaries = esummary_data.get("DocumentSummarySet", {}).get("DocumentSummary", [])
            if summaries:
                docsum = summaries[0]

        genomic_info = docsum.get("genomicinfo") or docsum.get("GenomicInfo") or []
        if isinstance(genomic_info, dict):
            genomic_info = [genomic_info]

        log_gene_lookup(f"docsum keys={list(docsum.keys()) if docsum else []}")

        if not genomic_info:
            log_gene_lookup("no genomicinfo in docsum")
            return None, "未在NCBI记录中找到基因坐标"

        region = genomic_info[0]
        chrom = (
            region.get("ChrLoc")
            or region.get("chr")
            or docsum.get("chromosome")
            or docsum.get("Chromosome")
        )
        start = (
            region.get("ChrStart")
            or region.get("chrstart")
            or docsum.get("chrstart")
        )
        end = (
            region.get("ChrStop")
            or region.get("chrstop")
            or docsum.get("chrstop")
        )

        log_gene_lookup(
            f"region keys={list(region.keys())}; raw chrom={chrom} start={start} end={end}"
        )

        if start is None or end is None or chrom is None:
            log_gene_lookup(
                f"missing fields after fallback chrom={chrom} start={start} end={end}; docsum keys={list(docsum.keys())}"
            )
            return None, "NCBI返回数据不完整，缺少染色体或坐标"

        # 标准化染色体格式
        chrom_str = str(chrom)
        if chrom_str.upper() in ["MT", "M"]:
            chrom_str = "chrM"
        elif not chrom_str.lower().startswith("chr"):
            chrom_str = f"chr{chrom_str}"

        try:
            start_int = int(start)
            end_int = int(end)
        except Exception:
            start_int = int(float(start)) if start is not None else None
            end_int = int(float(end)) if end is not None else None

        start_pos = int(min(start_int, end_int))
        end_pos = int(max(start_int, end_int))
        center_pos = int((start_pos + end_pos) / 2)

        log_gene_lookup(f"parsed chrom={chrom_str} start={start_pos} end={end_pos} center={center_pos}")

        return {
            "query": gene_symbol,
            "gene_id": gene_id,
            "chromosome": chrom_str,
            "start": start_pos,
            "end": end_pos,
            "center": center_pos,
            "strand": region.get("ChrStrand"),
            "map_location": docsum.get("maplocation") or docsum.get("MapLocation"),
            "summary": docsum.get("summary") or docsum.get("Summary"),
        }, "查询成功"
    except requests.RequestException as req_err:
        log_gene_lookup(f"request error: {req_err}")
        return None, f"网络请求失败: {req_err}"
    except Exception as e:
        log_gene_lookup(f"parse error: {e}")
        return None, f"解析NCBI返回数据失败: {e}"


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fetch gene coordinates from NCBI and output TSV/JSON.")
    p.add_argument("-g", "--gene", action="append", help="Gene symbol, can be used multiple times")
    p.add_argument("-f", "--file", help="File with one gene symbol per line")
    p.add_argument("-o", "--organism", default="Homo sapiens", help="Organism name (default: 'Homo sapiens')")
    p.add_argument("--format", choices=["tsv", "json"], default="tsv", help="输出格式，tsv 或 json (默认 tsv)")
    return p.parse_args(argv)


def read_genes_from_file(path: str) -> List[str]:
    with open(path, "r") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    genes: List[str] = []
    if args.gene:
        for g in args.gene:
            genes.extend([s.strip() for s in g.split(",") if s.strip()])
    if args.file:
        genes.extend(read_genes_from_file(args.file))

    # 如果没有通过参数给出，从 stdin 读取（支持管道）
    if not genes and not sys.stdin.isatty():
        genes.extend([line.strip() for line in sys.stdin if line.strip()])

    if not genes:
        log_gene_lookup("No gene provided. Use -g/--gene, -f/--file, or pipe names into stdin.")
        print("请通过 -g 或 -f 提供基因，或通过管道输入基因列表", file=sys.stderr)
        return 2

    results = []
    for gene in genes:
        res, msg = fetch_gene_coordinates(gene, organism=args.organism)
        if res is None:
            results.append({"query": gene, "error": msg})
        else:
            results.append(res)

    if args.format == "json":
        print(json.dumps(results, ensure_ascii=False, indent=2))
    else:
        # TSV 输出
        fieldnames = ["query", "gene_id", "chromosome", "start", "end", "center", "strand", "map_location", "summary"]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())