# cNMF-Workthrough

### After starting up environment can navigate to the appropriate working directory by first backing out of the initial directory twice:  <br/> cd .. <br/> cd .. <br/> Can then go to the appropriate working directory with the following code: <br/> cd ./mnt/c/Programming_Stuff/cNMF-Workthrough ### 

Can download the data by entering the following url:<br/>
https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

Packages installed in Mamba environment.

Python: 
Version 3.10.16;
mamba install -c conda-forge python=3.10.16

Scikit-learn:
Version 1.3.2;
mamba install -c conda-forge scikit-learn=1.3

scanpy:
Version 1.9.8;
mamba install -c bioconda scanpy=1.9

anndata:
This was already installed with a previous package as version 11.3 which should work according to the GitHub for the cNMF package. 

R:
Version 4.4.2 "Pile of Leaves"
mamba install -c conda-forge r-base=4.4.2
### Originally tried installing 4.3.3 "Angel Food Cake" which was automatically updated to 4.4.2 "Pile of Leaves" after installing Seurat v5.1.0 ###

Seurat:
Version 5.1.0;
mamba install -c conda-forge r-seurat="5.1.0"

cnmf:
Version 1.6.0;
mamba install bioconda::cnmf

tidyverse
Version 2.0.0;
mamba install conda-forge::r-tidyverse

openxlsx
Version 4.2.8;
mamba install conda-forge::r-openxlsx


When unzipping the .gz.tar file using the system command from R I kept getting an error regarding not being able to find the file. I think this is due to using Windows subsystem for Linux. Was able to unzip, uncompress,and remove the gz.tar file with the following code <br/> gunzip filename.tar.gz <br/> tar -xf filename.tar <br/> rm filename.tar
 
