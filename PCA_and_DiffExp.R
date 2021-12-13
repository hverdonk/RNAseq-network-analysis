# load required packages
library(dplyr)
library(Seurat)
library(RColorBrewer)
library(SeuratObject)
library(SeuratWrappers)
library(DESeq2)
library(plotly)
library(BiocParallel)
register(MulticoreParam(detectCores(-1)))

# Enter directory name here
path_to_dir <- ""
setwd(path_to_dir)

# Read count data (use read.csv for csv file)
count_file_name <- ""
counts <- read.table(count_file_name)

# Read metadata (sample names, embryonic stage, sex)
metadata_file_name <- ""
meta <- read.csv(metadata_file_name)
meta$sex_stage <- as.factor(paste(meta$sex, meta$stage, sep="_"))

# Clean up count data 
names(counts) <- seq(0,ncol(counts)-1,1)

# Create Seurat object to filter low reads
s <- CreateSeuratObject(counts, min.cells = 3, min.features = 350)

# Create DESeq2 object for variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(s@assays$RNA@counts, 
                              colData = meta, 
                              design = ~sex_stage)

# Apply variance stabilizing transformation 
vst <- vst(dds, blind = FALSE)

# PCA 
pca <- plotPCA(vst, intgroup=c("sex","stage"))

######## Differential expression analysis #########
# Run Differential Expression 
dds <- DESeq(dds)

# Find differentially expressed genes between group1 and group2
# (group1 and group2 should be male and female groups from the same embryonic stage)
group1 <- ""
group2 <- ""
contrast=c("sex_stage", group1, group2)

res <- results(dds, contrast = contrast)

# Shrink imprecise fold changes
res <- lfcShrink(dds, contrast = contrast, res = res, 
                 type = "normal", parallel = TRUE)

# Export results 
output_filename <- ""
write.table(res, output_filename)
