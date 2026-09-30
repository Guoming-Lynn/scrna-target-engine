# Changelog

## 0.1.0 - 2026-09-30

- 从速览工作区拷出引擎代码，不包含 `hadha_source`。
- 富集改为只读本地 GMT。没有在线富集调用。
- 全局聚类报告写明分辨率由启发式自动选择。
- `BATCH_COLUMN` 的默认值只在 `config_contract.py`，三个阶段都从那里读。
- `inspect_obs.py` 只读 h5ad 的 obs 列名和少量取值。
- 本地 ORA 用调用方传入的基因背景和 scipy 超几何。GSEA 改为对本地 GMT 调用 `gseapy.prerank`，同一种子可重复。代码里没有 `.enrichr(`。
- 拆开超过 400 行的 phase01 和 plot_utils。卫生检查和 Windows/Linux CI 跑 pytest。`legacy/` 不进测试收集，也不进 wheel。
- Scrublet 使用配置里的 `RANDOM_SEED`。ORA 的 BH 分母包含与背景有交集但重叠为 0 的基因集。
