import os
import argparse
import pysam
from pymongo import MongoClient, ASCENDING
from tqdm import tqdm
import time

def get_mongo_client(host='127.0.0.1', port=30001):
    """建立数据库连接"""
    try:
        # 显式使用 127.0.0.1 避免 localhost 解析问题
        client = MongoClient(f'mongodb://{host}:{port}/', serverSelectionTimeoutMS=5000)
        client.server_info() # 触发一次连接检查
        return client
    except Exception as e:
        print(f"[Error] 无法连接到 MongoDB ({host}:{port})。请检查容器是否运行。")
        print(f"详细错误: {e}")
        exit(1)

def ingest_bam(bam_file, sample_name, db_name="worf_db", host='127.0.0.1', port=30001):
    # 1. 连接数据库
    client = get_mongo_client(host, port)
    db = client[db_name]
    collection = db['reads_mapping']

    print(f"[*] 正在处理样本: {sample_name}")
    print(f"[*] 读取 BAM 文件: {bam_file}")
    
    # 2. 读取 BAM
    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
    except ValueError:
        print("[Error] BAM 文件索引可能损坏或无法读取。")
        return

    # 3. 准备批量写入
    BATCH_SIZE = 50000
    batch_data = []
    total_inserted = 0
    
    # 仅仅为了进度条估算（可选，如果太慢可以去掉 total）
    # total_reads = samfile.count() 
    # print(f"[*] 预计处理 Reads: {total_reads}")
    
    start_time = time.time()
    
    # 4. 核心循环
    # 如果安装了 tqdm，这里会显示进度条
    for read in tqdm(samfile.fetch(), desc="Ingesting", unit=" reads"):
        if read.is_unmapped:
            continue
            
        # 修正染色体名称 (去掉 chr 前缀或保持一致，根据你的习惯调整)
        # 这里保持原样
        chrom = read.reference_name
        
        # 构建文档
        doc = {
            "s": sample_name,          # s: sample_id (简写省空间)
            "r": read.query_name,      # r: read_name
            "c": chrom,                # c: chrom
            "st": read.reference_start, # st: start (0-based)
            "ed": read.reference_end,   # ed: end (0-based, exclusive)
            "rl": read.query_length,    # rl: read_length (bases)
            "ml": read.reference_length,# ml: mapped_length (on genome)
            "mq": read.mapping_quality, # mq: mapq
            # "cig": read.cigarstring   # 可选：如果不需要展示CIGAR，注释掉以节省空间
        }
        
        batch_data.append(doc)
        
        # 达到批次大小，写入
        if len(batch_data) >= BATCH_SIZE:
            collection.insert_many(batch_data)
            total_inserted += len(batch_data)
            batch_data = []

    # 写入剩余数据
    if batch_data:
        collection.insert_many(batch_data)
        total_inserted += len(batch_data)
        
    samfile.close()
    
    print(f"\n[Success] {total_inserted} 条 Reads 已入库。耗时: {time.time()-start_time:.2f}s")
    
    # 5. 建立索引 (最关键的一步)
    # 只有在第一次入库或显式要求时才检查/建立索引，防止重复操作拖慢速度
    # 这里我们检查一下是否已有索引，没有则建立
    existing_indexes = collection.index_information()
    
    # 索引 A: 位置查询核心索引 (样本 -> 染色体 -> 起始 -> 终止)
    if "loc_idx" not in existing_indexes:
        print("[*] 正在构建位置索引 (这可能需要几分钟)...")
        collection.create_index(
            [("s", ASCENDING), ("c", ASCENDING), ("st", ASCENDING), ("ed", ASCENDING)],
            name="loc_idx",
            background=True # 后台建立，不阻塞数据库
        )
        print("    -> 位置索引构建指令已发送。")
        
    # 索引 B: Read ID 查询索引
    if "read_idx" not in existing_indexes:
        print("[*] 正在构建 Read ID 索引...")
        collection.create_index(
            [("r", ASCENDING)],
            name="read_idx",
            background=True
        )
        print("    -> Read ID 索引构建指令已发送。")

    client.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="将 BAM 文件 Reads 信息导入 MongoDB")
    parser.add_argument("--bam", required=True, help="BAM 文件路径")
    parser.add_argument("--sample", required=True, help="样本名称 (Sample ID)")
    parser.add_argument("--port", type=int, default=30001, help="MongoDB 端口 (Default: 30001)")
    parser.add_argument("--host", default="127.0.0.1", help="MongoDB Host IP")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.bam):
        print(f"[Error] 文件不存在: {args.bam}")
        exit(1)
        
    ingest_bam(args.bam, args.sample, port=args.port, host=args.host)