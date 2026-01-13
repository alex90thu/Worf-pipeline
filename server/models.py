from pydantic import BaseModel, Field
from typing import List, Optional

# 定义一条 Read 的返回格式
class ReadItem(BaseModel):
    read_name: str = Field(..., alias="r")
    chrom: str = Field(..., alias="c")
    start: int = Field(..., alias="st")
    end: int = Field(..., alias="ed")
    read_len: int = Field(..., alias="rl")
    map_len: int = Field(..., alias="ml")
    mapq: int = Field(..., alias="mq")
    # cigar: Optional[str] = Field(None, alias="cig") 

    class Config:
        populate_by_name = True  # 允许使用 alias (数据库里的 st) 填充

# 定义样品信息的返回格式
class SampleInfo(BaseModel):
    name: str
    ingest_time: str
    total_reads: int
    chrom_count: int