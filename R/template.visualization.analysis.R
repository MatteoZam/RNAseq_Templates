### Hands-on RNAseq - Visualization Template ###

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(EnhancedVolcano)

##### (1) Load Differential Expression Results and Metadata #####

# Remember to match paths and verify the row/column structures
res_file        <- "<file_name_e.g_results_df.rds_or_txt>"
vsd_file        <- "<file_name_e.g_vsd_obj.rds>"
annotation_file <- "<file_name_e.g_sample_ann.txt>"

# res_int_df <- readRDS(res_file)
# vsd        <- readRDS(vsd_file)
# sample_ann <- read.table(annotation_file, header = TRUE, sep = "\t")



##### (2) Define Publication-Ready Theme #####

# Custom ggplot2 theme tailored for scientific publishing layout
theme_science <- function(){ 
  theme_bw() + 
    theme(
      panel.grid = element_blank(),                                 # Remove background grid lines
      text = element_text(color = "black"),                         # Set default text color to black
      axis.text = element_text(color = "black"),                    # Set axis tick label color to black
      panel.border = element_rect(fill = NA, colour = "black", size = 1) # Enhance main bounding box outline
    )
}



##### (3) Volcano Plot #####

# Interactive/Customizable Volcano plots using EnhancedVolcano

# Define significance thresholds
threshold_p   <- 0.05
threshold_lfc <- 1

# Generate Volcano Plot object
volcano_plot <- EnhancedVolcano(res_int_df,
                                lab = rownames(res_int_df),
                                x = 'log2FoldChange',
                                y = 'padj',                    # Swap with 'pvalue' if padj is unavailable
                                title = 'Volcano Plot',
                                pCutoff = threshold_p,
                                FCcutoff = threshold_lfc,
                                pointSize = 3.0,
                                labSize = 6.0) +
                theme_science()
	
# Save volcano plot to output directory
pdf("volcano_plot_DEGs.pdf", width = 7, height = 7)
print(volcano_plot)
dev.off()



##### (4) Heatmap Data Preparation #####

# Isolate significant genes using absolute log2 fold-change and adjusted p-values
sig_genes <- rownames(res_int_df[which(res_int_df$padj < threshold_p & abs(res_int_df$log2FoldChange) > threshold_lfc), ])

# Extract stabilized counts matrix for selected features
mat <- assay(vsd)[sig_genes, ]

# Row-center and scale expression levels to compute Z-scores across samples
mat_scaled <- t(scale(t(mat)))



##### (5) Heatmap Metadata Annotations #####

# Subset specific column variables of interest from the sample metadata
selected_cols <- c("Neoplasm", "Gender", "Batch")
annotation_df <- sample_ann[, selected_cols, drop = FALSE]

# Ensure the annotation indices strictly equal the matrix column identifiers
rownames(annotation_df) <- sample_ann$ID

# Verify order alignment prior to rendering
identical(colnames(mat_scaled), rownames(annotation_df))

# Assign explicit hex/named palettes to discrete metadata classes
ann_colors <- list(
  Neoplasm = c("AML" = "red", "MDS" = "blue"),
  Gender   = c("Male" = "darkgreen", "Female" = "orange"),
  Batch    = c("Batch1" = "purple", "Batch2" = "cyan", "Batch3" = "grey")
)



##### (6) Heatmap Generation & Export #####

# Build custom diverging gradient (Blue = Downregulated, White = Unchanged, Red = Upregulated)
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Render and save the pheatmap plot directly 
pheatmap(mat_scaled, 
         color = my_palette,
         annotation_col = annotation_df,       # Map column annotation tracking
         annotation_colors = ann_colors,        # Inject defined cluster colors
         show_rownames = FALSE,                 # Hide cluttered row/gene labels
         show_colnames = FALSE,                 # Hide sample names for large sets
         cutree_rows = 6,                       # Split row dendrogram into 6 distinct blocks
         main = "Heatmap of DEGs: AML vs MDS",  # Plot header text
         filename = "heatmap_DEGs_AML_vs_MDS.pdf", # Target file path destination
         width = 8,                             # Physical document width constraints
         height = 10                            # Physical document height constraints
)