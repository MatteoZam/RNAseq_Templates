### hands-on RNAseq
library(omicsTools) ###
library(pcaMethods) ###

##### (1) Read count data and annotation file #####
# remeber to check the match between samples in counts and in sample annotation file
count_file <- "<file_name>"
sample_annotation_file <- "<file_name>"

# counts <- readRDS(count_file)
# counts <- read.table(count_file)
# sample_ann <- read.table(sample_annotation_file)

# check whether the order of samples in the dataset and annotation file is the same;
identical(colnames(counts), sample_ann$ID) #this assumes that ID is the column with sample IDs
if(!identical(colnames(counts), sample_ann$ID)){
  idx <- match(colnames(counts), sample_ann$ID)
  sample_ann <- sample_ann[idx, ]
}
identical(colnames(counts), sample_ann$ID)

##### (2) Library size #####
#Calculate library size and store it as a vector.
lib_size <- colSums(counts)
p <- plot_lib_size(counts) 
plot(p)

### (3)	FILTER low quality samples based on library size
# see the arguments subset.remove to remove samples and min.samples for the minimum numner of samples in which a gene must be detected
counts <- filter_counts(X = counts, filter.by.cpm = T, min.cpm = 3)

# if you remove some samples, update the annotation
identical(colnames(counts), sample_ann$ID) #this assumes that ID is the column with sample IDs
if(!identical(colnames(counts), sample_ann$ID)){
  idx <- match(colnames(counts), sample_ann$ID)
  sample_ann <- sample_ann[idx, ]
}
identical(colnames(counts), sample_ann$ID)


#### (4) Library size after filtering #####
#Calculate library size and store it as a vector.
lib_size <- colSums(counts)
p <- plot_lib_size(counts) 
plot(p)


#### (5) Normalization  #####
# if necessary, there are the optional arguments gene.annotation sample.annotation and design
norm_data_obj <- normalize_counts(X = counts)

#### (6) RLE ####
RLE_res <- RLE(SummarizedExperiment::assay(norm_data_obj$vsd))
RLE_res <- RLE(norm_data_obj$lcpm) #or on vst data using SummarizedExperiment::assay(norm_data_obj$vsd)

p <- plot_RLE(RLE_res)
plot(p)

#### (7) PCA ####
pca_res <- pca(t(SummarizedExperiment::assay(norm_data_obj$vsd)), nPcs = 5, scale = "uv")
pca_res <- pca(t(norm_data_obj$lcpm), nPcs = 5, scale = "uv")
