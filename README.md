# supervised Capacity Preserved Mapping (supCPM):  a supervised visualization method for scRNA-seq data
This repository provides the implementation of supCPM in our paper. supCPM is a supervised visualization methods, which could result a more interpretable figure.   

## Descriptions

**data**: contains raw data used in the paper.

**result**: contians results for the four examples in the paper.

**supCPM**: contains source codes for supCPM and run the supCPM.

**Plot for results.R**: plot figures in the paper.

**RunSeurat.R** : run Seurat on the four datasets for preprocess. 


## Tutorial
We introduce the basic idea of tuning the parameter and use an example in our paper, RNAmix, to illustrate how to use supCPM. 

### Parameters
**data**: a matrix of the input n*p dataset, where n is the cell number and p is the number of genes.

**label**: a vector indicates the label for each cell. 

**no_dims**: a scalar of the embedding dimensions. For visualization purpose, we choose 2 or 3. 

**compel_force**:  takes the value 0 or 1.  compel_force=1 if user wants to pull clusters a bit apart; compel_force =0 if user wants the best preservation of the geometry.

**geodist**: takes the value 0 or 1. 0 means the user wants to use Euclidean distance in the high diemensions, and 1 means the users wants to choose geodesic distance.  
According to our experience,  we recomend to choose the Euclidean distance for clean data (less noisy) and choose geodesic distance for the noisy one.    

**degree**: a scalar to choose the degree of freedom of high dimensions t-distribution. Normally, the value is set to be 1 if geodesic distance is chosen and 
is set to 2 if Euclidean distance is used. Due to the supervised part of supCPM, 
the cluster's shape tends to shrink. This parameter could in some degree prevents the cluster from streching into a long shape. 

**ratio**: takes the value from 0 to 1. A larger value means higher proportion of the supervised term. 

**k**: choose the size of knn neighbors used in geodesic distance. This parameter will not influence the results of supCPM with Euclidean distance.  The default value is set to 7.

**change**: the number of iteration for CPM part. 

**niter**: total number of iteration.  (If niter < change, there will be no supervision.)

**seed**: the random seed for reproduction. 

**factor**: the constant multiplies on the high dimensional iter-cluster distances. This parameter takes the value >1. However, it should not be too large, 
otherwise the geometry structure will be distorted. In all the four datasets, set it to 1.3 will be enough.  

### RNAmix
#### Preprocess
We use the workflow of Seurat (v4.0.1) on the RNAmix dataset including procedures like normalization, scaling, PCA, clustering and visualization. We follow the standard procedures of 
Seurat, which could also be found in their toutrial (\web)
```r
library(Seurat)
# loading the dataset, the variable 'dir' is the path of file 'data'. 
rnamix.data <- read.csv(paste0(dir,"/data/RNAmix/RNAmix_celseq2.count.csv"))
# Create the Seurat Object 
RNAmix.seurat <- CreateSeuratObject(counts = rnamix.data)
# Normalize the dataset
RNAmix.seurat <- NormalizeData(RNAmix.seurat)
RNAmix.seurat <- FindVariableFeatures(RNAmix.seurat, selection.method = "vst", nfeatures = 2000)
all.genes  <- rownames(RNAmix.seurat)
# Scale the dataset
RNAmix.seurat <- ScaleData(RNAmix.seurat, features = all.genes)
# Run PCA
RNAmix.seurat <- RunPCA(RNAmix.seurat, features = VariableFeatures(object = RNAmix.seurat))
# Run SNN based clustering
RNAmix.seurat <- FindNeighbors(RNAmix.seurat, dims = 1:10)
RNAmix.seurat <- FindClusters(RNAmix.seurat, resolution = 1)
# Run UMAP
RNAmix.seurat <- RunUMAP(RNAmix.seurat, dims = 1:10)
# Run t-SNE
RNAmix.seurat <- RunTSNE(RNAmix.seurat, dims = 1:10)
# these procedures are also aggregated into a function called RunSeurat in our RunSeurat.R file.
# the code 'RNAmix.seurat <- RunSeurat(rnamix.data, npcs = 10, norm='log', resolution = 1, nfeature=2000)' yeilds the same results.
```

#### supCPM and CPM
With the labels from SNN based clustering in Seurat and first 10 PCs, we are ready to run the code of supCPM and CPM.  The implementation of other visualization methods used in our paper could be found 
in the file 'RunCPM.m' in detail. We use supUMAP from <https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap> 
and supPCA from Prof. Ali Ghodsi. 

Run the CPM (niter > change, factor=1 and degree =1)
``` 
rnamix_cpm = supCPM(rnamix_pca(:,1:10),rnamix_label,2,1,0,1,0,8,1000,900,123,1);
```
Run the supCPM
```
rnamix_supCPM_eu = supCPM(rnamix_pca(:,1:10),rnamix_label,2,1,0,2,0.7,7,500,2000,123,1.3);
rnamix_supCPM_geo = supCPM(rnamix_pca(:,1:10),rnamix_label,2,1,1,1,0.7,7,400,2000,123,1.3);
```

