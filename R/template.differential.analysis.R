library(omicsTools) ###

# REQUIREMENTS
# sample_ann: data.frame with sample annotation
# dge: DGE object; it can be obtained using normalize_counts()
# design: design matrix that can be obtained by means of model.matrix()
out_dir <- "<output_directory>"

### LIMMA
cont.matrix <- limma::makeContrasts(HIGHvsLOW=DISEASE_CATEGORYHR_MDS - DISEASE_CATEGORYLR_MDS, levels=design)
degs_limma <- differential_expression_limma(dge = dge, design = design, contr_mat = cont.matrix, out_dir = out_dir)


### DESeq
degs_DESeq <- differential_expression_DESeq(counts = counts, sample_ann = sample_ann, design = ~ DISEASE_CATEGORY + Sex + Run, contrast = c("DISEASE_CATEGORY", "AML", "CTRL"))


