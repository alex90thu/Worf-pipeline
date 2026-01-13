from pymongo import MongoClient
from datetime import datetime

# 配置
MONGO_URL = "mongodb://127.0.0.1:30001/"
DB_NAME = "worf_db"

class DBHandler:
    def __init__(self):
        self.client = MongoClient(MONGO_URL)
        self.db = self.client[DB_NAME]
        self.col = self.db['reads_mapping']

    def get_samples(self):
        pipeline = [
            {"$group": {"_id": "$s", "first_id": {"$first": "$_id"}}},
            {"$sort": {"_id": 1}}
        ]
        results = list(self.col.aggregate(pipeline))
        samples = []
        for res in results:
            t = res['first_id'].generation_time.astimezone()
            # 简单统计一下该样品的总数
            count = self.col.count_documents({"s": res['_id']})
            samples.append({
                "name": res['_id'],
                "ingest_time": t.strftime("%Y-%m-%d %H:%M:%S"),
                "total_reads": count,
                "chrom_count": 0 # 暂时不查详细的，为了API响应速度
            })
        return samples

    def query_reads(self, sample: str, chrom: str, start: int = None, end: int = None, limit: int = 100):
        query = {"s": sample, "c": chrom}
        if end: query["st"] = {"$lt": end}
        if start: query["ed"] = {"$gt": start}
        
        # 返回迭代器，并在转换时映射字段
        cursor = self.col.find(query).sort("st", 1).limit(limit)
        return list(cursor)

# 创建一个全局实例
db = DBHandler()