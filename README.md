# cNMF-Workthrough

### After starting up environment can navigate to the appropriate working directory by first backing out of the initial directory twice:  <br/> cd .. <br/> cd .. <br/> Can then go to the appropriate working directory with the following code: <br/> cd ./mnt/c/Programming_Stuff/cNMF-Workthrough ### 


Packages installed in Mamba environment.

Python: 
Version 3.10.16;
mamba install -c conda-forge python=3.10.16

Scikit-learn:
Version 1.4.2;
mamba install -c conda-forge scikit-learn=1.4

scanpy:
Version 1.8.2;
mamba install -c bioconda scanpy=1.8

anndata:
Version 0.11.3;
mamba install conda-forge::anndata

R:
Version 4.4.2 "Pile of Leaves"
mamba install -c conda-forge r-base=4.4.2
### Originally tried installing 4.3.3 "Angel Food Cake" which was automatically updated to 4.4.2 "Pile of Leaves" after installing Seurat v5.1.0 ###

Seurat:
Version 5.1.0;
mamba install -c conda-forge r-seurat="5.1.0"
 
