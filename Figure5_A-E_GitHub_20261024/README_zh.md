# Figure 5A–E GitHub 数据包

本包根据 2026-10-24 提供的最终 Figure 5 编号整理。原目录 43 份数据 CSV 全部原样复制，文件内容和精度不变；未重新运行差异表达或富集分析。

## 当前编号

| 当前图 | 上传目录 | 原目录 |
|---|---|---|
| 5A 火山图 | data/A_volcano | A_volcano |
| 5B 全部差异基因热图 | data/B_DEG_heatmap | B_DEG_heatmap |
| 5C 功能程序热图 | data/C_functional_program_heatmap | E_functional_program_heatmap |
| 5D 成脂相关转录因子 | data/D_adipogenic_TF_heatmap | F_Adipogenic_TF_heatmap |
| 5E KEGG | data/E_KEGG | C_KEGG |

## 上传哪些文件

将本目录内的文件及子目录整体作为 GitHub 仓库内容，保留目录结构。ZIP 是本地传输副本，解压后使用即可，无需在仓库里再重复放一份 ZIP。

- `data/shared`：完整输入 counts、原样本信息、size factors。
- `data/A_volcano` 至 `data/E_KEGG`：每个图的表达数据、统计结果及实际绘图数值。
- `supporting_analysis`：同一已完成分析中的完整受检基因表达矩阵、GSEA 全结果与选取依据、过滤/映射审计，补足图中筛选结果的背景。
- `README.md`、`DATA_DICTIONARY.md`、`metadata`：英文使用说明、字段解释、原分析配置和软件环境。
- `provenance`、`validation`、`manifest.csv`：来源、旧编号对照和核对证据。
- `code`：只读核验脚本及依赖；不是新的差异表达分析。

## 核对结论

43/43 份数据 CSV 的原 SHA-256 校验通过；原目录核对 155 项全部通过。A 的 17,882 个点与 10 个标签，B 的 897 上调/926 下调基因，C 的 24 个基因及 FDR，D 的 14 个基因及 9 个星号，E 的上下调各 10 条通路均与最终图内容相符。详细范围见 `AUDIT_REPORT_zh.md`。

完整输入 counts 有 77,775 行；17,887 行进入模型。5 个 NA 检验结果未画在火山图中，已保留在完整表及补充矩阵。原始 counts 指计数矩阵，本包不含 FASTQ/BAM；不能把 VST、Z 分数或标准化 counts 称为原始测序数据。

原组名 NC 在当前图中显示为 Vehicle。所有数据列名保留原样，通过 `metadata/sample_display_mapping.csv` 说明。`0*` 是下溢脚注；D 基因名后的 `*` 表示全基因 BH-FDR <0.05。

## 中文核对：公开发布时补充

补充作者/联系方式、GitHub 仓库网址、论文标识（如已有）及作者确定的数据和代码许可。本地数据整理已完成，尚未上传 GitHub；未擅自指定许可或生成 accession/DOI。作者若另有 FASTQ 或测序数据库记录，应在 README 中关联。
