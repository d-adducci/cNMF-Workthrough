
#==============================================================================#
# Following through the example code from the cNMF github for running the 
# package through R. This allows augmentation of Seurat analyses. 
# A lot of the comments in this file are directly from the tutorial. 
# www.github.com/dylkot/cNMF/blob/master/Tutorials/R_vignette.Rmd
#==============================================================================#

library(data.table) # Version 1.16.4
library(Matrix) # Version 1.7.2
library(Seurat) # Version 5.1.0
library(tidyverse) # Version 2.0.0
library(openxlsx) # Version 4.2.8

# Where data is stored
data_dir <- "./ProcessedData/"

# Where the filtered data will be stored 
filtered_dir <- "./ProcessedData/filtered/"

# Name of the run being performed
runname <- "example_cNMF"

# Where the tar.gz file is located
#file_path <- paste0(data_dir,"pbmc_unsorted_3k_filtered_feature_bc_matrix.tar.gz")


#==============================================================================#
# Data Setup
#==============================================================================#

# I am assuming 'system' command pushes the code written in R to the 
# environment where it can be ran in Python/non-R 
# This code ended up not working. Pretty sure it is do to running on Windows
# Was able to unzip and uncompress file manually from the shell. 

# Extracting tar.gz file
#extract_command <- sprintf("tar -xzf %s -C %s", shQuote(file_path), 
#                           shQuote(raw_dir))
#system(extract_command) 

# Removing the zip file after extracting
# Do not need this for the way I did it
#rm_command <- sprintf("rm %s", file_path)
#system(rm_command)

# Storing data locally
#system("ls /filetered_gene_bc_matrices/hg19", intern = TRUE)

# Data can then be loaded in with the following code
pbmc.data <- Read10X(data.dir="./ProcessedData/filtered_gene_bc_matrices/hg19")



#==============================================================================#
# Loading in and formatting unzipped data
#==============================================================================#

# Loading in raw data
pbmc.data <- Read10X(data.dir="./ProcessedData/filtered_gene_bc_matrices/hg19")

# Initializing the Seurat object with the raw (non-normalized data) 
# Some filtering performed which removed genes detected in less than 3 cells
# and cells with less than 200 gene counts.
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nCount_RNA > 200)

# Checking dimensions of Seurat object
dim(pbmc)

# Extracting necessary information from Seurat object
counts <- pbmc@assays$RNA$counts
barcodes <- colnames(counts)
gene_names <- rownames(counts)

# Checking counts and barcodes
counts[1:5, 1:5]
barcodes[1:5]

#------------------------------------------------------------------------------#
# Saving the filtered and extracted information to the disk
#------------------------------------------------------------------------------#

filtered_dir <- "./ProcessedData/filtered/"

if(!dir.exists(filtered_dir)){
  dir.create(filtered_dir, recursive = TRUE)
}

# Output counts matrix
writeMM(counts, paste0(filtered_dir,"matrix.mtx"))

# Output cell barcodes
barcodes <- colnames(counts)
write.table(as.data.frame(barcodes),paste0(filtered_dir,'barcodes.tsv'),
            col.names = FALSE, row.names = FALSE, sep = "\t")

# Output feature names
gene_names <- rownames(counts)
features <- data.frame("gene_id" = gene_names,"gene_name" = gene_names,
                       type = "Gene Expression")
write.table(as.data.frame(features),sep = "\t",paste0(filtered_dir,'genes.tsv'),
            col.names = FALSE, row.names = FALSE)

#==============================================================================3
# Running cNMF
#==============================================================================#

# First step is the prepare step which normalizes the count matrix and 
# prepares the factorization step
runname <- "example_cNMF"
cmd <- paste("cnmf prepare --output-dir", data_dir,
             "--name", runname,
             "-c", paste0(filtered_dir,'matrix.mtx'),
             "--max-nmf-iter 2000",
             "-k 5 6 7 8 9 10 --n-iter 20", sep=" ")
print(cmd)
system(cmd)


# Next is the factorization sep which runs the NMF --n-iter times (which 
# is set to 20 in this case) for each value of K. All jobs are run 
# sequentially here on a single worker, however these can be distributed to
# multiple cores or nodes using separate commands such as below

# cnmf factorize --output-dir ./R_Example_Data/ --name example_cNMF --worker-index 0 --total-workers 3
# cnmf factorize --output-dir ./R_Example_Data/ --name example_cNMF --worker-index 1 --total-workers 3
# cnmf factorize --output-dir ./R_Example_Data/ --name example_cNMF --worker-index 2 --total-workers 3

# Factorization
cmd <- paste("cnmf factorize --output-dir", data_dir,
            "--name", runname,
            "--worker-index 0 --total-workers 1",sep = " ")
print(cmd)
system(cmd)


# Next we concatenate the results for each value of K into a single file
cmd <- paste("cnmf combine --output-dir", data_dir,
             "--name", runname, sep = " ")
print(cmd)
system(cmd)


# A plot is made estimating the trade-off between higher value of K 
# and stability and error
cmd <- paste("cnmf k_selection_plot --output-dir", data_dir,
             "--name", runname, sep = " ")
print(cmd)
system(cmd)


# The plot suggest K=7 might be a local optima in stability 
# (maximize stability, minimize error) So that will be the first option explored
cmd <- paste("cnmf consensus --output-dir", data_dir,
             "--name", runname,
             "--components", 7,
             "--local-density-threshold", 0.1,
             "--show-clustering", sep = " ")
print(cmd)
system(cmd)

#==============================================================================#
# Adding to Seurat to identify GEPs
#==============================================================================#


# Consensus scoring looks clean, indicating that 7 is a good starting point
usage_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "usages", "k_7.dt_0_1","consensus","txt",sep="."),sep="/")
spectra_score_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "gene_spectra_score", "k_7.dt_0_1","txt",sep="."),sep="/")
spectra_tpm_file <- paste(data_dir[1:length(data_dir)], runname, paste(runname, "gene_spectra_tpm", "k_7.dt_0_1","txt",sep="."),sep="/")

usage <- read.table(usage_file, sep="\t", row.names=1, header=TRUE)
spectra_score <- read.table(spectra_score_file, sep="\t", row.names=1, header=TRUE)
spectra_tpm <- read.table(spectra_tpm_file, sep="\t",row.names=1, header=TRUE)
head(usage)


# For most analyses the resulting usage file output is normalized so 
# that each cell sums to 1. Code below is what does that
usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))


# Concatenating the usage_norm into the metadata of the Seurat object
# to make it easier to plot later
new_metadata <- merge(pbmc@meta.data, usage_norm, 
                      by = "row.names", all.x = TRUE)
rownames(new_metadata) <- new_metadata$Row.names
pbmc@meta.data <- new_metadata


# Running the standard Seurat UMAP pipeline so the GEP usages can be 
# plotted over the UMAP
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = all.genes)
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- RunUMAP(pbmc, dims = 1:15)


# Setting up directory to save images
image_dir <- "./Plots/"

# Plotting the learned GEPs
p <- FeaturePlot(pbmc, features = colnames(usage_norm), combine = F)

pdf(file=paste0(image_dir,"UMAP_GEPs.pdf"))
print(p)
dev.off()


# Extracting the top 20 most highly weighted genes for each GEP below
get_top_colnames <- function(row){
  # Orders the values in descending order and gets the names of the top 20
  print(row[1:5])
  top_indices <- order(row, decreasing = TRUE)[1:20]
  return(colnames(spectra_score)[top_indices])
}

top_colnames <- apply(spectra_score, 1, get_top_colnames)
top_colnames <- as.data.frame(top_colnames)

write.xlsx(x=top_colnames,file="./ProcessedData/top_20_spectra_scores.xlsx")




