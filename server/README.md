# Worf Genomics API 说明

本服务基于 FastAPI 提供对 MongoDB 中 WGS 映射数据的查询接口。

- 服务名称：Worf Pipeline Server
- 框架：FastAPI
- 数据库：MongoDB（默认连接 `mongodb://127.0.0.1:30001/`）
- 数据库名：`worf_db`
- 集合：`reads_mapping`

## 启动方式

确保已安装依赖（FastAPI、Uvicorn、pymongo 等）。在项目根目录运行，根据自己的情况调整端口：

```bash
uvicorn server.main:app --reload --host 0.0.0.0 --port 8018
```

启动后可访问：`http://localhost:8018`

交互式文档（Swagger UI）：`http://localhost:8018/docs`

## 认证与权限

当前接口未启用认证，默认公开。若对外提供服务，建议在网关层或应用层增加认证与访问控制。

## 接口概览

### 1) 健康检查
- 方法：GET
- 路径：`/`
- 说明：服务在线状态检查
- 响应示例：
```json
{
  "status": "online",
  "system": "Worf Pipeline Server"
}
```

### 2) 获取样品列表
- 方法：GET
- 路径：`/samples`
- 标签：Samples
- 说明：返回所有已入库样品的简要信息
- 响应模型：`SampleInfo[]`
- 字段：
  - `name`：样品名
  - `ingest_time`：入库时间（取集合中该样品第一条记录的 ObjectId 生成时间）
  - `total_reads`：该样品对应的记录总数
  - `chrom_count`：染色体数量（当前为占位，固定 0）
- 响应示例：
```json
[
  {
    "name": "SAMPLE_001",
    "ingest_time": "2025-12-31 23:59:59",
    "total_reads": 123456,
    "chrom_count": 0
  }
]
```

### 3) 查询 Reads（核心接口）
- 方法：GET
- 路径：`/query/{sample_id}`
- 标签：Query
- 说明：按样品与染色体位置范围查询映射到该区域的 Reads。
- 查询参数：
  - `chrom`（必填）：染色体名，如 `chr1`
  - `start`（选填，int，0-based）：起始位置（不含）
  - `end`（选填，int，exclusive）：结束位置（不含）
  - `limit`（选填，int，默认 100，范围 1~5000）：最大返回记录数
- 响应模型：`ReadItem[]`
- 字段说明（当前默认按 alias 输出短字段名）：
  - `r`：`read_name`（Reads 名称）
  - `c`：`chrom`（染色体）
  - `st`：`start`（起始坐标）
  - `ed`：`end`（结束坐标）
  - `rl`：`read_len`（Read 长度）
  - `ml`：`map_len`（Read 与参考的比对长度）
  - `mq`：`mapq`（比对质量分数）

> 说明：模型 `ReadItem` 使用了 Pydantic 的 `alias` 与数据库字段映射（如 `st`、`ed`）。FastAPI 的默认配置会按 alias 输出响应。如果希望返回人类可读的字段名（如 `read_name`、`start` 等），可在路由层设置 `response_model_by_alias=False`。

- 响应示例（按 alias 输出）：
```json
[
  {
    "r": "read_0001",
    "c": "chr1",
    "st": 12345,
    "ed": 12456,
    "rl": 111,
    "ml": 111,
    "mq": 60
  }
]
```

## 使用示例

- 获取样品列表：
```bash
curl -s "http://localhost:8018/samples" | jq .
```

- 查询 chr1 上指定区间的 Reads：
```bash
curl -s "http://localhost:8018/query/SAMPLE_001?chrom=chr1&start=100000&end=200000&limit=500" | jq .
```

- 查询整条染色体（受 `limit` 限制）：
```bash
curl -s "http://localhost:8018/query/SAMPLE_001?chrom=chr1&limit=1000" | jq .
```

### 命令行工具（worf_tools.sh）
- 如需在命令行完成入库或交互式查询，可使用根目录的 `worf_tools.sh`。
- 示例：
  - 打开交互式菜单：`./worf_tools.sh`
  - 手动入库：`./worf_tools.sh ingest <bam_path> [sample_id]`
  - 交互式查询（TUI）：`./worf_tools.sh query`
- 依赖：MongoDB (`127.0.0.1:30001`)，以及 `scripts/bam2mongo.py`、`scripts/query_worf.py`。

## 数据库结构与查询逻辑（简述）

- 集合：`reads_mapping`
- 关键字段（与 API 响应 alias 对应）：
  - `s`：样品名（对应 `sample_id`）
  - `c`：染色体
  - `st`：起始位置
  - `ed`：结束位置
  - 其他：`r`、`rl`、`ml`、`mq` 等
- 查询条件：
  - 必须匹配 `s=sample_id` 与 `c=chrom`
  - 若提供 `start`：匹配 `ed > start`
  - 若提供 `end`：匹配 `st < end`
  - 结果按 `st` 升序排序，并受 `limit` 限制

## 注意事项

- 若查询结果为空，返回空列表而非 404。
- `limit` 上限为 5000，谨慎提高以避免大响应体与数据库压力。
- 当前 `chrom_count` 未实时统计，后续可扩展为统计每个样品的染色体覆盖情况。

## 维护与扩展建议

- 若希望返回人类可读字段名，建议在 `server/main.py` 的路由上设置 `response_model_by_alias=False`。
- 如需增加过滤条件（如 `mapq` 阈值），可在 `server/database.py` 的 `query_reads()` 中扩展 `query` 字段，并在路由添加相应 Query 参数。
- 建议在生产环境中配置认证（Token/Key）和速率限制。
