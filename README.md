SupCPM:  supervised visualization method for scRNAseq data
---------------------------------------------------------------------------

# Introduction
supCPM (supervised Capacity Preserving Mapping) provides a clustering guided visualization algorithm, which could result a more interpretable figure. 

To find out the details of the supCPM, you could check out our preprint in BioRxiv,

[Zhiqian Zhai, Yu L. Lei, Rongrong Wang, Yuying Xie, **Supervised Capacity Preserving Mapping: A Clustering Guided Visualization Method for scRNAseq data**. *BioRxiv 2021.06.18.448900*](https://doi.org/10.1101/2021.06.18.448900)

And we greatly appreatiate Dr. Elnaz Barshan and Dr. Ali Ghodsi for providing us the codes of SPCA. [[Barshan, E. et al. (2011). Supervised principal component analysis: Visualization, classification and regression on subspaces and submanifolds. Pattern Recognition, 44(7), 1357–1371.](https://doi.org/10.1016/j.patcog.2010.12.015) ]

# Contents
* [Contents of Files in the Respository](#contents-of-files-in-the-respository)
* [R](#r)
  * [Installation](#installation)
  * [Functions](#functions)
  * [Parameters](#parameters)
  * [Quick Start](#quick-start)
  * [supCPM_downsample to accelerate the large dataset](#supcpm_downsample-to-accelerate-the-large-dataset)
  * [Tuning Parameters](#tuning-parameters)

* [MATLAB](#matlab)
  * [Paramters](#parameters-1)
  * [Example: RNAmix](#example-rnamix)
    * [Preprocess](#preprocess)
    * [Run supCPM and CPM](#run-supcpm-and-cpm)

# Contents of Files in the Repository
**data**:  Raw data we used in our paper.

**supCPM_MATLAB**: Contains implementation and results of supCPM in MATLAB.

**supCPM**: Contains the implementation of supCPM in R.

**supCPM figures**: Contains 1) main figures in our paper; 2) codes used in our paper including preprocess; running PCA, supCPM, t-SNE, UMAP; calculating merics and plotting figures.

**readme figures**: Contains figures used in this readme file.

# R 
## Installation
You could install supCPM directly from R with `devtools`:

    if (!suppressWarnings(require(devtools))) install.packages("devtools")
    devtools::install_github("zhiqianZ/supCPM/supCPM", depedencies=T)

If it doesn't work, you could also download the `.tar.gz` file from our github and manually install the package. In this way you may need to install the dependent packages manually as well (igraph, pracma, DescTools).

## Functions
Our package `supCPM` contains three main functions. The first is `supCPM`, which is the implementation of our algorithm. The second is `VisualMetric`, which contains five comparison metics used in our article. The third is `supCPM_downsample`. This function is to deal with the large sample dataset, which is more and more prevalent in scRNA, and it is still under development. Also, from the workflow of supCPM, it is easy to see that if we set some part of the parameters in supCPM, supCPM could result in CPM. (`degree`=1, `ratio`=any value, `niter2`=0, `factor`=1)
 
## Parameters
All the parameters used in the function `supCPM` are listed and explained in detail as follows:

`data` : The data matrix whose **rows are cells and columns are genes**

`label`: A label vector indicating which cluster the cell belongs to. The vector could be numeric with cluster ID or with characters identifying names of cell types. 

`ratio`: A numeric value in the objective function balancing off the geometric and label information. Default is 0.7.

`no_dims`: An integer suggest which dimensions the data needs to be projected into. Default value is 2.

`comel_force`: 0 or 1 to declare whether to pull clusters a bit apart during the intrinsic dimensional estimation. Default is 1. 

`dist`: A character to choose the distance measure in the high dimensions, 'Euclidean' or 'geodesic'. The parameter is **not case sensitive**. Default is 'Euclidean'. 

`degree`: A numeric value used to determine the degree of freedom of t-distribution in the high dimensions. As a generalization, this `degree` is not necessary to be an integer.  Default is 2.

`k`: An integer used to build a KNN graph to calculate geodesic distance. Thus, this parameter will **only influence** the result if you choose the `dist` to be 'geodesic'. Default is 7.

`niter1`:  An integer to determine the iteration runs for the first phase of supCPM. Default is 500.

`niter2`: An integer to determine the iteration runs for the second phase of supCPM. Default is 700.

`seed`: Random seed to generate the initialization of embedding points. If you choose `init`=True, then the initialization is fixed, and `seed` will not work. Default is 40.

`factor`: A numeric to multiple on the capacity-adjusted distance matrix to separate different clusters. Default is 1.3.

`verbose`: Whether to output the objective function value every ten interactions. Default is True.

`init`: Whether to use MDS coordinates of capacity adjusted distance as initialization. Default is False.

`epsilon`: The factor multiples on the default &epsilon<sub>0</sub> used in the high dimensional t-distribution to increase (>1) or decrease (<1) the original value. Default is 1.

`lr`: Learning rate. Default is 500.

`intermediate`: Whether to output the intermediate result or not. The default is False.

## Quick Start
Say if you have loaded a data matrix `RNAmix_pca`, storing the first 10 PCs of RNAmix (after performing normalization and scaling in Seurat) and a vector `RNAmix_label` (from SNN-based clustering in Seurat) in R. You can run supCPM with Euclidean distance and geodesic distance as follows:

    library(supCPM)
    RNAmix_supCPM <- supCPM(RNAmix_pca,RNAmix_label,ratio=0.6)
    RNAmix_supCPM <- supCPM(RNAimx_pca,RNAmix_label,raito=0.6,dist='geodesic',degree=1)
    
Then, if we plot the result, we could get a figure similar to this using `ggplot2`.

<img src="/readme figures/supCPM_demo.PNG" width="70%" />

## supCPM_downsample to accelerate the large dataset
We implement supCPM_downsample() to function as an alternative if your dataset has a large number of cells. 

It has three steps:

 (1) Sketching: use the [geometric sketching](https://github.com/brianhie/geosketch) to downsample a subset of cells from a large dataset.
 
 (2) supCPM: supCPM() on the downsampled cells
 
 (3) Projection: project the remaining cells onto the downsampled supCPM space as the centers of their kNN among the subsampled cells plus a random perturbation

Compared with supCPM(), supCPM_downsample() has two additional parameters: (1) `alpha` (default as 0.1) and (2) `return.idx` (default as FALSE). 

`alpha` represents the percentage of cells that are subsampled to actually run the supCPM algorithm; 

`return.idx` indicates whether to return the cell-wise label showing whether the cell is subsampled to run supCPM (labeled as 'sketched cell') or projected to the subsampled supCPM space (labeled as 'projected cell'). If `return.idx` = TRUE, it will return a list of two elements: supCPM coordinates named 'supCPM' and downsampled labels named 'idx'.

**Prior to running supCPM_downsample(), you need to install the Python package 'sketch' and the R package 'reticulate.'** 
Here is a short demo for the installation and running supCPM_downsample:

    install.packages("reticulate")
    library(reticulate)
    library(supCPM)
    conda_create("supCPM-sketch")
    use_condaenv("supCPM-sketch")
    py_install("geosketch", envname="supCPM-sketch",pip=T)
    supcpm = supCPM_downsample(pca, label, alpha=0.1, niter1 = 300, niter2 = 300, ratio=0.5, return.idx=T)
    supCPM = supcpm$supCPM
    skectch_label = supcpm$idx

## Tuning parameters
There are several parameters that need to be tuned in supCPM. Luckily, most of these parameters can be judged and chosen by looking at the visualization results. The details and examples can be found in the supplementary materials of our publication. Here, we discuss how to choose them.

`ratio`: The ratio is the most important parameter for supCPM to supervise. It works to pull cells within the same cluster to be closer to each other and repel cells from a different cluster. So larger ratio means you want clusters to be more dense and scattered. We recommend that users try from 0.6-0.7 and increase it if there are still incorrect placements. This could be checked by plotting cells with labels. Also, we argue that the smallest ratio that could separate clusters is the best way to avoid further distortion. 

`factor`: The factor coarsely categorizes clusters in the high dimensions. By multiplying a factor on the capacity-adjusted distance matrix, clusters are separated slightly. Notice that 1) both `factor` and `ratio` are needed. `factor` operates overly, while there might still be a small fraction of cells residing in wrong clusters due to the noise of scRNA data. That is where `ratio` takes effect. Setting `factor` to be a large value would separate different clusters clearly. It seems `ratio` is not necessary. However, increasing inter-cluster distances without limitation is a problem not only for supCPM but also for UMAP and t-SNE. Extremely large distances will be converted into a small probability with minor distinction. So, the algorithms could fail to preserve the long distances. Thus, if we can't set `factor` to be large, then we need to refine cell positions individually by `ratio`. Empirically, `factor`= 1.3 is good enough.  The following figures shows a brief example of why we need both `ratio` and `factor`.
<img src="/readme figures/ratio&factor_demo.PNG" width="100%" />

`degree`: The degree of freedom in the t-distribution plays a role in preventing cluster shape from shrinking anisotropically. KL-divergence, as a part of the objective function, will impose different forces in different directions based on the distribution of points. So when the algorithm tries to shrink the cluster, it will do anisotropically.  This will result in clusters with long shapes. So, we need to decrease the inter-cluster influence by increasing the probability. Set `degree` to 2 for Euclidean distance and 1 for geodesic distance would suit most cases. If you witness the result of plenty long shape clsuters, which doesn't happen before the switching of obejctive function, this suggests you to increase the parameter `degree`. The effect of the degree of freedom can be seen in the following figures. This is a toy example, where three spheres on a line are visualized in supCPM with different degrees of freedom.
<img src="/readme figures/degree_demo.PNG" width="100%" />

`epsilon`: The `epsilon` is used in the high dimensional t-distribution to prevent the inverse of zero. However, the magnitude of it could influence the structure of data. If the `epsilon` is too small, the cluster will lose its original shape, and the density of the points will look artificial. These phenomena indicate us to tune `\epsilon` a little bit larger. If `\epsilon` is too large, then the probability similarity will be smoothed. Thus, clusters are squeezed, and only a relatively large distance can be captured. Observing this pattern means the necessity of decreasing `\epsilon`.
<img src="/readme figures/epsilon.PNG" width="100%" />

#  MATLAB
We have introduced the basic idea of tuning the parameters in the previous chapter, and here we just simply summarize the parameters in the MATLAB version, which are slightly different parameters. Then, use an example in our paper, RNAmix, to illustrate how to use supCPM. 
## Parameters
`data`: a matrix of the input n*p dataset, where n is the cell number, and p is the number of genes.

`label`: a vector indicates the label for each cell. 

`no_dims`: a scalar of the embedding dimensions. For visualization purposes, we choose 2 or 3. 

`compel_force`: takes the value 0 or 1. compel_force=1 if the user wants to pull clusters a bit apart; compel_force =0 if the user wants the best preservation of the geometry.

`geodesist`: takes the value 0 or 1. 0 means the user wants to use Euclidean distance in the high dimensions, and 1 means the user wants to choose the geodesic distance.  
Based on our experience, we recommend choosing the Euclidean distance for clean data (less noisy) and the geodesic distance for noisy data.    

`degree`: a scalar to choose the degree of freedom of high dimensions t-distribution. Normally, the value is set to be 1 if geodesic distance is chosen and 
is set to 2 if Euclidean distance is used. Due to the supervised part of supCPM, 
the cluster's shape tends to shrink. This parameter could, to some degree, prevent the cluster from stretching into a long shape. 

`ratio`: takes the value from 0 to 1. A larger value means a higher proportion of the supervised term. 

`k`: choose the size of known neighbors used in the geodesic distance. This parameter will not influence the results of supCPM with Euclidean distance.  The default value is set to 7.

`change`: the number of iterations for the CPM part. 

`niter`: total number of iterations.  (If niter < change, there will be no supervision.)

`seed`: the random seed for reproduction. 

`factor`: the constant multiplies on the high dimensional iter-cluster distances. This parameter takes the value >1. However, it should not be too large.
Otherwise the geometry structure will be distorted. In all four datasets, setting it to 1.3 will be enough.  

# Example: RNAmix
## Preprocess
We use Seurat's workflow on the RNAmix dataset, including procedures like normalization, scaling, PCA, clustering, and visualization. We follow the standard procedures of 
Seurat, which could also be found in their tutorial of PBMC3k <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>.
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

## Run supCPM and CPM
With the labels from SNN-based clustering in Seurat and the first 10 PCs, we are ready to run the code of supCPM and CPM. The settings for other visualization methods used in our paper can be found in the file 'RunCPM.m' in detail. We use supUMAP from <https://www.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap> 
and the codes of supPCA are provided by SPCA authors Dr. Ali Ghodsi and Dr. Elnaz Barshan. 

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


