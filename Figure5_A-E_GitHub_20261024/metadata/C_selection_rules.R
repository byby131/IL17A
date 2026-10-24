# Rules fixed for this user-requested, post-discovery Figure 5 extension.
# This is NOT pre-registration and does NOT alter the original inferential families.
figure5_rules <- list(
  global_programs_each_direction=5L,
  programs_each_function=2L,
  leading_edge_genes_each_program=6L,
  program_jaccard_threshold=.5,
  program_order="FDR_global_HR ascending, abs(NES) descending, set_id ascending",
  gene_order="abs(Wald statistic) descending, abs(log2FC) descending, detectable samples descending, gene_id ascending",
  focus_patterns=c(
    adipocyte_differentiation="ADIPOGENESIS|ADIPOCYTE_DIFFERENTIATION",
    lipid_storage_turnover="TRIGLYCERIDE|TRIACYLGLYCER|LIPID_DROPLET|LIPID_STORAGE|LIPOLYSIS|LIPOGENESIS",
    lipid_metabolic_processes="FATTY_ACID|FATTY_ACYL|LIPID_METABOLISM|CHOLESTEROL|PLASMA_LIPOPROTEIN|PPAR",
    mature_function_related_not_specific="INSULIN_RECEPTOR|GLUCOSE_TRANSPORT|GLUT4|ADIPONECTIN|LEPTIN"
  ),
  function_patterns=c(
    lipid_function="ADIPOGENESIS|ADIPOCYTE_DIFFERENTIATION|TRIGLYCERIDE|TRIACYLGLYCER|LIPID_DROPLET|LIPID_STORAGE|LIPOLYSIS|LIPOGENESIS|FATTY_ACID|FATTY_ACYL|LIPID_METABOLISM|CHOLESTEROL|PLASMA_LIPOPROTEIN|PPAR",
    inflammation_chemokines="TNFA_SIGNALING_VIA_NFKB|^HALLMARK_INFLAMMATORY_RESPONSE$|^REACTOME_CHEMOKINE|^REACTOME_CYTOKINE",
    ECM_organization="EXTRACELLULAR_MATRIX|COLLAGEN|ECM"
  ),
  explicit_existing_readouts="PPARG",
  TGF_response_annotation="Hallmark::HALLMARK_TGF_BETA_SIGNALING",
  TGF_BMP_search="TGF|BMP|BONE_MORPHOGEN|SMAD|TRANSFORMING_GROWTH",
  early_data_status="transcript_counts_received_quantification_and_design_metadata_pending",
  font_family="Arial",font_face="bold",font_pt=10,png_dpi=600
)
