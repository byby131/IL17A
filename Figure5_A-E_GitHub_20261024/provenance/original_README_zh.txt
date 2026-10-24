Figure5_A4_editable：五个面板的原始counts及绘图源数据

对应最终文件 Figure5_A4_editable.pdf（已确认v9版）。本次仅拆分已有数据，未重新拟合或改变统计结果。
所有CSV可用Excel、R或Python打开。新拆分CSV为UTF-8 BOM。原始表格副本保留原文件内容。

A_volcano：17,882个实际绘制散点、全部17,887个受检基因统计、原图10个标注基因，以及绘制基因的counts/标准化表达/VST。
B_DEG_heatmap：全部1,823个DEG，897上调、926下调；按图中从上到下顺序保存。
C_KEGG：上调、下调各10条通路，附横条长度、命中基因、命中基因counts及统计、富集背景和全量富集结果。
E_functional_program_heatmap：原图全部24个功能程序基因，含行注释和程序颜色对应表。
F_Adipogenic_TF_heatmap：原图全部14个成脂相关转录因子。
shared：最初输入的完整counts文件、样本分组及DESeq2 size factors。

表达矩阵六列均按NC_1、NC_2、NC_3、IL17_1、IL17_2、IL17_3排列。
NC为第7天MDI+vehicle；IL17为MDI+IL-17A 50 ng/mL，各3个独立生物学重复。
raw_counts是原始输入counts；normalized_counts是DESeq2标准化counts；VST_expression是方差稳定变换后的表达，不能作为原始counts使用。
heatmap_VST_rowZ是实际热图使用的逐基因Z分数。heatmap_color_values_clipped仅反映颜色饱和范围：B为[-2.5,2.5]，E/F为[-2,2]，未截断数值另表完整保留。
row_order为热图从上到下顺序；draw_order为火山图散点绘制顺序。
log2FC为IL17/NC的未收缩log2倍数变化；正值表示IL17组高。baseMean为平均标准化counts；lfcSE为log2FC标准误；stat为Wald统计量；pvalue为原始P值；padj为BH校正FDR。
DEG定义：padj<0.05，无额外倍数变化阈值。NA是缺失值，不是0。源数据FDR=0属于数值下溢，对应图中0*，没有替换成新的P值。
A的plot_x_log2FC截断于[-8,8]；plot_y_minus_log10_FDR截断于60。minus_log10_FDR在源FDR为0时按绘图代码使用float最小正规正数作为计算下限。这些均为显示变换，原始统计列保持不变。
C的bar_length_minus_log10_FDR为横条长度；Count为命中基因数；GeneRatio=Count/foreground_size；BgRatio=pathway_size/universe_size；FDR为原富集结果的多重检验校正值。上下调是在各自DEG列表中富集，不等价于证明通路激活/抑制。
KEGG的hit_genes和hit_entrez分别是集合列表，不应按字符串位置一一配对。pathway_gene_membership以基因符号展开，保留一基因参与多条通路的关系。
E/F图中已删除Log2FC显示列，但源数据中的效应量完整保留。row_annotations保留源分析注释，其中未在最终PDF显示的注释不作为该图额外结论。
程序P1-P5分别为Fatty acid metabolism、TNF/NF-kB signaling、Chemokine receptors、ECM organization、Collagen formation。完整定义见E文件夹。F图星号表示FDR<0.05。

manifest.csv记录每份CSV的来源、行数和SHA256校验值。CSV直接打开Excel可能触发基因名称自动日期转换；建议通过“数据-从文本/CSV”导入并将gene_id列指定为文本。
