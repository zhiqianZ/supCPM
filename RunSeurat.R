dir = '...' # where the file locate 
result.dir = '/result/output/'
library(Seurat)
library(Matrix)
library(tidyverse)
# self-defined functions
RunSeurat <- function(count,npcs,norm='sct',resolution=0.8,nfeature=3000){
  if(norm=='sct'){
    seurat.obj <- CreateSeuratObject(counts = count)
    seurat.obj <- PercentageFeatureSet(seurat.obj, pattern = "^MT-", col.name = "percent.mt")
    seurat.obj <- SCTransform(seurat.obj, vars.to.regress = "percent.mt", verbose = FALSE,
                              variable.features.n = nfeature)
    seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npcs, verbose=F)
    seurat.obj <- FindClusters(seurat.obj, resolution = resolution, verbose = F)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npcs,verbose=F)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npcs,verbose=F)
  }
  if(norm=='log'){
    seurat.obj <- CreateSeuratObject(counts = count)
    seurat.obj <- NormalizeData(seurat.obj)
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = nfeature)
    all.genes  <- rownames(seurat.obj)
    seurat.obj <- ScaleData(seurat.obj, features = all.genes)
    seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
    seurat.obj <- FindNeighbors(seurat.obj, dims = 1:npcs)
    seurat.obj <- FindClusters(seurat.obj, resolution = resolution)
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:npcs)
    seurat.obj <- RunTSNE(seurat.obj, dims = 1:npcs)
  }
  return(seurat.obj)
}
save.result <- function(seurat.obj,dir,name,norm='sct'){
  write.csv(seurat.obj@reductions$tsne@cell.embeddings,paste0(dir,name,'_tsne.csv'))
  write.csv(seurat.obj@reductions$umap@cell.embeddings,paste0(dir,name,'_umap.csv'))
  write.csv(seurat.obj@reductions$pca@cell.embeddings,paste0(dir,name,'_pca.csv'))
  if(norm=='sct'){
    write.csv(seurat.obj@assays$SCT@scale.data,paste0(dir,name,'_scale.csv'))
  }else{
    write.csv(seurat.obj@assays$RNA@scale.data,paste0(dir,name,'_scale.csv'))
  }
}
load_mixseq_data <- function(expt_dir,condition) {
  counts <- Matrix::readMM(file.path(expt_dir, 'matrix.mtx'))
  rownames(counts) <- readr::read_tsv(file.path(expt_dir, 'genes.tsv'), col_names = F)$X2
  counts <- counts[!duplicated(rownames(counts)),]
  colnames(counts) <- readr::read_tsv(file.path(expt_dir, 'barcodes.tsv'), col_names = F)$X1
  classifications <- readr::read_csv(file.path(expt_dir, 'classifications.csv')) %>%
    dplyr::mutate(condition = condition) %>% tibble::column_to_rownames("barcode") %>%
    as.data.frame()
  sc_data <- Seurat::CreateSeuratObject(counts,meta.data = classifications)
}
###################### Synthetic #######################
nn <- c(600, 600, 600, 400, 400, 400, 200, 200, 200)
d1 <- c(5,   5,   5,   8,  8,  8,   10,  10,  10)
d2 <- c(20,  20,  20,  20,  20,  20,  20,  20, 20)
cl <- c(0,   0,   0,   1,   1,   1,   2,   2,   2)
v <- c(1,1,1, 1.2,1.2,1.2, 2,2,2)
p <- 50
X <- NULL
synthetic_label <- NULL
set.seed(123)
for (i in 1:9){
  Xpart = v[i] * matrix(rnorm(nn[i]*p),ncol=50)
  Xpart[,i] = Xpart[,i] + d1[i]
  Xpart[,20+cl[i]] = Xpart[,20+cl[i]] + d2[i]
  X <- rbind(X,Xpart)
  ypart = rep(i,nn[i])
  synthetic_label <- c(synthetic_label,ypart)
}
rownames(X) <- paste('point',1:nrow(X))
colnames(X) <- paste('dim',1:ncol(X))
synthetic.mean <- colMeans(X)
X.centered <- t(X) - synthetic.mean

# Create SeuratObject
synthetic <- CreateSeuratObject(counts = t(X))
synthetic <- FindVariableFeatures(synthetic, selection.method = "disp",nfeatures = 50)
synthetic@assays$RNA@scale.data <- X.centered
synthetic <- RunPCA(synthetic,features = VariableFeatures(object = synthetic))
synthetic <- RunUMAP(synthetic,dims=1:20)
synthetic <- RunTSNE(synthetic,dims=1:20)



###################### PBMC3K #######################
pbmc.3k <- Read10X(data.dir= paste0(dir,"/data/pbmc3k/filtered_gene_bc_matrices/hg19/"))
pbmc <- RunSeurat(pbmc.3k,npcs=30)
# Annotate Celltypes 
new.cluster.ids <- c("Naive CD4 T","Memory CD4 T","CD14+ Mono","B","FCGR3A+Mono","Effector CD8 T",
                     "NK","Memory CD8 T","Naive CD8 T","DC","IFN-activated CD4T","Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

###################### RNAmix #######################
rnamix.data <- read.csv(paste0(dir,"/data/RNAmix/RNAmix_celseq2.count.csv"))
rnamix <- RunSeurat(rnamix.data, npcs = 10, norm='log', resolution = 1, nfeature=2000)


##################### Cancer #######################
cancer_data <- load_mixseq_data(".../data/rna/DMSO_24hr_expt1","cancer")
# filter cells
cancer_meta <- cancer_data@meta.data
cell_quality <- cancer_meta['cell_quality']
normal <- which(cell_quality=='normal')
cancer_data <- cancer_data@assays$RNA@counts
cancer_data_filtered <- cancer_data[,normal]
cancer_label <- cancer_meta['singlet_ID']
cancer_label_filtered <- cancer_label[normal,]
cancer_label_filtered <- cbind(cancer_label_filtered,as.numeric(as.factor(cancer_label_filtered)))
# select clusters
cancer_part <- cancer_data_filtered[,cancer_label_filtered[,2]%in%c(11,12,17,14,16,22,6,9,1,2,15)] 
cancer_label_part <- cancer_label_filtered[cancer_label_filtered[,2]%in%c(11,12,17,14,16,22,6,9,1,2,15),] 
# SCTransform
cancer <- RunSeurat(cancer_part,npcs = 20,resolution = 0.4)
# rename the SNN clustering result (SCT)
celltype <- names(summary(as.factor(cancer_label_part[,1])))
idents <- celltype[c(1,8,5,11,10,3,7,9,4,2,6)]
names(idents) <- levels(cancer)
cancer <- RenameIdents(cancer, idents)
# accuracy 
correct = 0
sum(as.character(cancer@active.ident)==cancer_label_part[,1])/length(cancer_label_part[,1])

################ save the result ################
save.result(pbmc,paste0(dir,result.dir),"pbmc3k")
save.result(rnamix,paste0(dir,result.dir),"RNAmix",'log')
save.result(cancer,paste0(dir,result.dir),'cancer')

label.pbmc <- as.data.frame(pbmc@active.ident)
label.pbmc <- cbind(label.pbmc,as.numeric(pbmc@active.ident))
colnames(label.pbmc) <- c('celltype','cluster id')
label.cancer <- cbind(as.character(cancer@active.ident),as.numeric(cancer@active.ident))
write.csv(label.pbmc,paste0(dir,result.dir,'pbmc3k','_label.csv'))
write.csv(rnamix@active.ident,paste0(dir,result.dir,'RNAmix','_label.csv'))
write.csv(label.cancer,paste0(dir,result.dir,'cancer','_label.csv'))


write.csv(synthetic_label,paste0(dir,result.dir,'synthetic_label.csv'))
write.csv(synthetic@reductions$pca@cell.embeddings,paste0(dir,result.dir,'synthetic_pca.csv'))
write.csv(synthetic@reductions$tsne@cell.embeddings,paste0(dir,result.dir,'synthetic_tsne.csv'))
write.csv(synthetic@reductions$umap@cell.embeddings,paste0(dir,result.dir,'synthetic_umap.csv'))
write.csv(X,paste0(dir,result.dir,'synthetic.csv'))
