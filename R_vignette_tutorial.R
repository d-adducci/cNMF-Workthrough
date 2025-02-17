
#==============================================================================#
# Following through the example code from the cNMF github for running the 
# package through R. This allows augmentation of Seurat analyses. 
# A lot of the comments in this file are directly from the tutorial. 
# www.github.com/dylkot/cNMF/blob/master/Tutorials/R_vignette.Rmd
#==============================================================================#

library(data.table) # Version 1.16.4
library(Matrix) # Version 1.7.2
library(Seurat) # Version 5.1.0

# Where raw data is stored
data_dir <- "/RawData/"

# Where the tar.gz file is located
file_path <- paste0(data_dir,"pbmc_unsorted_3k_filtered_feature_bc_matrix.tar.gz")

#==============================================================================#
# Data Setup
#==============================================================================#

# I am assuming 'system' command pushes the code written in R to the 
# environment where it can be ran in Python/non-R 

# Extracting tar.gz file
extract_command <- sprintf("tar -xzf %s -C %s", shQuote(file_path), 
                           shQuote(data_dir))
system(extract_command) 

# Removing the zip file after extracting
rm_command <- sprintf("rm %s", file_path)
system(rm_command)

# Storing the data locally
system("ls ./RawData/pbmc_unsorted_3k_filtered_feature_bc_matrices/hg19",
       intern = TRUE)

# Data can then be loaded in with the following code
pbmc.data <- Read10X(data.dir="/RawData/pbmc_unsorted_3k_filtered_feature_bc_matrices/hg19",
                     intern = TRUE)



#==============================================================================#
# Loading in and formatting unzipped data
#==============================================================================#

# Loading in raw data
pbmc.data <- Read10X(data.dir="/RawData/pbmc_unsorted_3k_filtered_feature_bc_matrices/hg19",
                     intern = TRUE)

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

filtered_dir <- "/ProcessedData/"

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






