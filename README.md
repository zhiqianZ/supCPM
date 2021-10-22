supCPM:  supervised visualization method for scRNAseq data
---------------------------------------------------------------------------

## Introduction
supCPM (supervised Capacity Preserving Mapping) provides a clustering guided visualization algorithm, which could result a more interpretable figure. 

To find out the details of the supCPM, you could check out our preprint in BioRxiv,

[Zhiqian Zhai, Yu L. Lei, Rongrong Wang, Yuying Xie, **Supervised Capacity Preserving Mapping: A Clustering Guided Visualization Method for scRNAseq data**.]
(bioRxiv 2021.06.18.448900; doi: https://doi.org/10.1101/2021.06.18.448900)

### Contents
* [Contents of Files in the Respository] (#contents-of-files-in-the-respository)

* [R] (#r)
**[Installation](#installation)
**[Functions](#functions)
**[Parameters](#parameters-1)
**[Quick Start](#quick-start)
**[Tuning Parameters](#tuning parameters)

* [Matlab] (#matlab)
** [Paramters] (#parameters)
** [Example](#example)
*** [Preprocess](#preprocess)
*** [Run supCPM and CPM](#run-supcpm-and-cpm)

### Contents of Files in the Respository
**data**  Raw data we used in our paper.
**MATLAB** Contains implementation and results of supCPM in MATLAB.
**R** Contains implementaion of supCPM in R.

### R 
#### Installation
First, you could install supCPM directly from R with `devtools`:

    if (!suppressWarnings(require(devtools))) install.packages("devtools")
    devtools::install_github("zhiqianZ/supCPM")

If it doesn't work, you could also download the `.tar.gz` file from our github and manually install the package.

#### Functions
Our package `supCPM` contains three main functions. The first is `supCPM`, which is the implementation of our algorithm. The second is `VisualMetric`, which contains five comparison metics used in our article. The third is `supCPM_downsample`. This function is to deal with the large sample dataset, which is more and more prevelent in scRNA, and it is still under development. Also, from the workflow of supCPM, it is easy to see that if we set some part of the parameters in supCPM, supCPM could result in CPM. (`degree`=1, `ratio`=any value, `niter2`=0, `factor`=1)
 
#### Parameters
All the parameters used in the function `supCPM` is listed and explained in detail as following:

`data` : The data matrix whose **rows are cells and columns are genes**

`label`: A label vector indicating which cluster the cell belongs to. The vector could be numeric with cluster ID or with characters identifying names of cell types. 

`ratio`: A numeric value in the objective function balancing off the geometric and label information.

`no_dims`: An integer suggest which dimensions the data needs to be projected into. Default value is 2.

`comel_force**: 0 or 1 to declare whether to pull clusters a bit apart during the intrinsic dimensional estimation. Default is 1. 

`dist`: A character to choose the distance measure in the high dimensions, 'Euclidean' or 'geodesic'. The parameter is **not case sensitive**. Default is 'Euclidean'. 

`degree`: A numeric value to determine the degree of freedom of t-distribution used in the high dimensions. As a generalizaiton, this `degree` is not neccessary to be an integer.  Default is 2.

`k`: An integer used to build a KNN graph for the calculation of geodesic distance. Thus, this parameter will **only influence** the result if you choose the `dist` to be 'geodesic'. Default is 7.

`niter1`:  An integer to determine the interation runs for the first phase of supCPM. Default is 500.

`niter2`: An integer to determine the interation runs for the second phase of supCPM. Default is 700.

`seed`: Random seed to generate the initialization of embedding points. If you choose `init`=True, then the initialzaiton is fixed and `seed` will not work. Default is 40.

`factor`: A numeric to multiple on the capacity adjusted distance matrix to separate different clusters. Default is 1.3.

`verbose`: Whether to ouput the objective funciton value every ten interations. Also, the scatterplot for the embedding points before switching of objective function will be shown if `verbose`=True. Default is True.

`init`: Whether to use MDS coordinates of capacity adjusted distance as initialization. Default is False.

`epsilon`: The factor multiples on the default $\epsilon_0$ used in the high dimensional t-distribution to increase (>1) or decrease (<1) the original value. Default is 1.

`lr`: Learning rate. Default is 500.

`inter`: Whether to output intermediate result or not. Default is False. This could be used to tune the parameters as will be explained in the Tuning Parameters section.  

#### Quick Start
Say, if you have loaded a data matrix `RNAmix_pca`, storing the first a few PCs of RNAmix and a vector `RNAmix_label` in R. You can run supCPM with Euclidean distance and geodesic distance as follows:

    library(supCPM)
    RNAmix_supCPM <- supCPM(RNAmix_pca,RNAmix_label,ratio=0.6)
    RNAmix_supCPM <- supCPM(RNAimx_pca,RNAmix_label,raito=0.6,dist='geodesic',degree=1)

#### Tuning parameters
There are several parameters needs to be tuned in supCPM. But luckily, most of these parameters could be judged and chosen by looking at the visualization results. The details and examples could be found in the supplemnetary materials of our publication. Here, we simply discuss how to choose them.

`ratio`: The ratio is the most important parameter for supCPM to do the supervision.  Ratio works for pull cells within the same cluster to be closer with each other and repel cells from a different cluster. So larger ratio means you want clusters to be more dense and scattered. We recommand users to try from 0.6-0.7 and increase it if there are still incorrect placement. This could be check by plotting cells with labels. Also, we argue that the smallest ratio that could separate clusters is the best to aviod further distortion. 

`factor`: The factor is to separate clusters in the high dimensions coarsely. By multiplying a factor on the capacity adjusted distance matrix, clusters are separated a little bit. Notice that 1) both `factor` and `ratio` are needed. `factor` operates overly, while there might still be a small fraction of cells resides in wrong clusters due to the noise of scRNA data. That is where `ratio` takes effect. Setting `factor` to be a large value would separate different clusters clearly. It seems `ratio` is not necessary. But increasing inter-cluster distances without limitation is a problem for not only supCPM, but UMAP and t-SNE as well. Extremely large distances will be converted into small probability with minor distinction. So the algorithms could fail to preserve the long distances. Thus, if we can't set `factor` to be large, then we need to refine cell positions individually by `ratio`. Emperically, `factor`= 1.3 is good enough.   

`degree`: The degree of freedom in the t-distribution plays a role in preventing cluster shape from shrinking anisotrophically. KL-divergence, as a part of objective function, will impose different force in different direction based on the distribution of points. So when the algorithm try to shrink cluster, it will do anisotropically.  This will result in clusters with long shape. So what we need to do is to decrease the inter-cluster influence by increasing the probability. Set `degree` to 2 for Euclidean distance and to 1 for geodesic distance would be suitable for most cases. If you witness the result of plenty long shape clsuters, which doesn't happen before the switching of obejctive function, this suggests you to increase the parameter `degree`.

`epsilon`: The `epsilon` is used in the high dimensional t-distribution to prevent the inverse of zero. However, the magnitude of it could influence the structure of data overly. In short, if `epsilon` is too small, you will see the strange pattern that points on the boundary of cluster identically distributed on the circle. This tells you to increase `epsilon`. If `epsilon` is too large, you will see clusters in long shape before the switching of objective function, which is different from `degree` case. Then, you need to decrease `epsilon`. What's more, by setting `inter`=T, you could output both figure and matrix for the embedding cells before switching the objective function.


###  MATLAB
We introduce the basic idea of tuning the parameter and use an example in our paper, RNAmix, to illustrate how to use supCPM. 
#### Parameters
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

### Example RNAmix
#### Preprocess
We use the workflow of Seurat (v4.0.1) on the RNAmix dataset including procedures like normalization, scaling, PCA, clustering and visualization. We follow the standard procedures of 
Seurat, which could also be found in their toutrial of PBMC3k <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>.
```r
library(Seurat)
# Load the dataset where the variable 'dir' is the path of file 'data'. 
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
With the labels from SNN based clustering in Seurat and first 10 PCs, we are ready to run the code of supCPM and CPM.  The settings for other visualization methods used in our paper could be found 
in the file 'RunCPM.m' in detail. We use supUMAP from <https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap> 
and supPCA from Prof. Ali Ghodsi. 

Run the CPM (niter > change, factor=1 and degree =1)
``` 
rnamix_cpm = supCPM(rnamix_pca(:,1:10),rnamix_label,2,1,0,1,0,8,1000,900,123,1);
```
Run the supCPM
```
# Euclidean distance
rnamix_supCPM_eu = supCPM(rnamix_pca(:,1:10),rnamix_label,2,1,0,2,0.7,7,500,2000,123,1.3);
# Geodesic distance
rnamix_supCPM_geo = supCPM(rnamix_pca(:,1:10),rnamix_label,2,1,1,1,0.7,7,400,2000,123,1.3);
```


