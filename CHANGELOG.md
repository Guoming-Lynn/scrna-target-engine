# Changelog

## 0.1.0 - 2026-09-30

- 根目录许可证改为专有声明：Copyright 2026 GuomingLin. All rights reserved.
- 从速览工作区拷出引擎代码，不包含 `hadha_source`。
- 富集改为只读本地 GMT。没有在线富集调用。
- 全局聚类报告写明分辨率由启发式自动选择。
- `BATCH_COLUMN` 的默认值只在 `config_contract.py`，三个阶段都从那里读。
- `inspect_obs.py` 只读 h5ad 的 obs 列名和少量取值。
- 本地 ORA 用调用方传入的基因背景和 scipy 超几何。GSEA 改为对本地 GMT 调用 `gseapy.prerank`，同一种子可重复。代码里没有 `.enrichr(`。
- 拆开超过 400 行的 phase01 和 plot_utils。卫生检查和 Windows/Linux CI 跑 pytest。`legacy/` 不进测试收集，也不进 wheel。
- Scrublet 使用配置里的 `RANDOM_SEED`。ORA 的 BH 分母包含与背景有交集但重叠为 0 的基因集。
- `write_demo_h5ad.py` 写出一份带固定种子的合成矩阵，供 Studio 演示。不是策展数据。
- 合成矩阵可以在断网时跑完四个阶段。80 个细胞的图质量仍会因为小提琴面板细胞数不足被拒绝。空的经典 marker 点图和没有显著基因的差异表不再中断流程。Harmony 不再向当前 harmonypy 请求 1 个 cluster。
- 10x h5 用 `scanpy.read_10x_h5` 读入。测试用合成的两细胞矩阵，不使用真实 Cell Ranger 文件。
- 汇总图清单记下每一格的位置，以及旁边同名的 csv 或 tsv。
- 演示矩阵默认写成 1600 个细胞、250 个基因。80 个细胞仍可用 `--cells 80` 做快速回归，那次图质量会拒绝。
