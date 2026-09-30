# scrna-target-engine

细胞水平速览引擎。不做 donor 水平推断。

## 安装

Python 3.11。当前 pin 的 anndata 0.12.19 不支持 3.10。在仓库根目录：

```text
uv sync
```

依赖版本写在 `pyproject.toml`。富集读取本地 GMT。GSEA 使用 `gseapy.prerank`，不调用在线接口。

## 启动

准备一个本地 GMT 文件，在配置里写 `LOCAL_GMT`。没有这个文件时，差异基因仍会计算，富集记录为未运行，基因列表不会送出本机。

```text
python scripts/run_pipeline.py --config config.json
```

`legacy/phase04_pseudotime.py` 不在这条主流程里。

## 测试

```text
uv run pytest
```

分辨率报告包含「分辨率由启发式自动选择」。
