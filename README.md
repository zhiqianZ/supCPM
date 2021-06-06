# supervised Capacity Preserved Mapping (supCPM):  a supervised visualization method for scRNA-seq data
This repository provides the implementation of supCPM in our paper. supCPM is a supervised visualization methods, which could result a more interpretable figure.   

## Descriptions

**data**: contains raw data used in the paper.

**result**: contians results for the four examples in the paper.

**supCPM**: contains source codes for supCPM and run the supCPM.

**Plot for results.R**: plot figures in the paper.

**RunSeurat.R** : run Seurat on the four datasets for preprocess. 


## Tutorial
We use an example in our paper, RNAmix, to illustrate how to use supCPM  and introduce the basic idea of tuning the parameter. 

####Preprocess
supCPM requires the scaled dataset and labels as the input. 


####Parameters
**data**: a matrix of the input n*p dataset, where n is the cell number and p is the number of genes.

**label**: a vector indicates the label for each cell. 

**no_dims**: a scalar of the embedding dimensions. For visualization purpose, we choose 2 or 3. 

**compel_force**:  takes the value 0 or 1.  compel_force=1 if user wants to pull clusters a bit apart; compel_force =0 if user wants the best preservation of the geometry.

**geodist**: takes the value 0 or 1. 0 means the user wants to use Euclidean distance in the high diemensions, and 1 means the users wants to choose geodesic distance.  
According to our experience,  we recomend to choose the Euclidean distance for clean dataet (less noisy) and choose geodesic distance for noisy dataset.    

**degree**: a scalar to choose the degree of freedom of high dimensions t-distribution. Normally, the value is set to be 1 if geodesic distance is chosen and 
is set to 2 if Euclidean distance is used. Due to the supervised part of supCPM, 
the cluster's shape tends to shrink. This parameter could in some degree prevents the cluster from streching into a long shape. 

**ratio**: takes the value from 0 to 1. A larger value means higher proportion of the supervised term. 

**k**: choose the size of knn neighbors used in geodesic distance. This parameter will not influence the results of supCPM with Euclidean distance. 

**change**: the number of iteration for CPM part. 

**niter**: total number of iteration.  (If niter < change, there will be no supervision.)

**seed**: the random seed for reproduction. 

**factor**: the constant multiplies on the high dimensional iter-cluster distances. This parameter takes the value >1. However, it should not be too large, 
otherwise the geometry structure will be distorted. In all the four datasets, set it to 1.3 will be enough.  
