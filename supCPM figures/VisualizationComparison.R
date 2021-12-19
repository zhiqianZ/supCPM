dir = '.../SupCPM-main/' #where the raw date locates
result.dir = '.../supCPM-main/supCPM_MATLAB/results/'  # results for supUMAP, MDS and supPCA needs to load from matlab results
fig.dir = "..."  # where your figure want to save
library(cowplot)
library(Seurat)
library(Matrix)
library(tidyverse)
library(supCPM)
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
synthetic.label <- NULL
set.seed(123)
for (i in 1:9){
  Xpart = v[i] * matrix(rnorm(nn[i]*p),ncol=50)
  Xpart[,i] = Xpart[,i] + d1[i]
  Xpart[,20+cl[i]] = Xpart[,20+cl[i]] + d2[i]
  X <- rbind(X,Xpart)
  ypart = rep(i,nn[i])
  synthetic.label <- c(synthetic.label,ypart)
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
synthetic.pca      <- synthetic@reductions$pca@cell.embeddings
synthetic.umap     <- synthetic@reductions$umap@cell.embeddings
synthetic.tsne     <- synthetic@reductions$tsne@cell.embeddings
synthetic.cpm      <- supCPM(synthetic.pca[,1:20],synthetic.label,dist='geodesic',degree=1,factor=1,niter2=0,init=F)
synthetic.supcpm_euclidean  <- supCPM(synthetic.pca[,1:20],synthetic.label,ratio=0.6,init=T)
synthetic.supcpm_geodesic <- supCPM(synthetic.pca[,1:20],synthetic.label,ratio=0.7,dist='geodesic',degree=1,init=T)
synthetic.supumap <- read.csv(paste0(result.dir,'synthetic_supUMAP.csv'),header=F)
synthetic.suppca <- read.csv(paste0(result.dir,'synthetic_supPCA.csv'),header=F)
synthetic.mds <- read.csv(paste0(result.dir,'synthetic_MDS.csv'),header=F)
## metrics
plot.title <- c("MDS","PCA","supPCA","tSNE","UMAP","supUMAP","CPM","supCPM_Euclidean","supCPM_Geodesic")
metric.data.frame <- data.frame(matrix(rep(0,160*4),ncol=4))
colnames(metric.data.frame) <- c('dataset','method','metric','Scores')
k<-1
for(i in 2:9){
   metric = eval(str2lang(paste0("synthetic.",tolower(plot.title[i]))))%>%VisualMetric(synthetic.pca[,1:20],synthetic.label,2,10)
   for(j in 1:5){
   metric.data.frame[k,1:3] = c('Synthetic',plot.title[[i]],names(metric)[[j]])
   metric.data.frame[k,4] = metric[[j]]
   k<-k+1
   }
}



## Plot
color.synthetic <- rgb(matrix(c(0,136,64,99,189,95,167,218,86,255,74,70,
                                163,0,89,255,46,127,27,229,255,0,110,165,
                                0,0,165),ncol=3,byrow = T),maxColorValue = 255)
synthetic.result <- list()
for(i in 1:9){
  result <- eval(str2lang(paste0("synthetic.",tolower(plot.title[i]))))
  colnames(result) <- paste0('V',1:ncol(result))
  synthetic.result[[i]] <- as.data.frame(result)
}
# Plot for synthetic data
# set seed for jittering
set.seed(123)
synthetic.plot <- list()
len <- max(abs(diff(range(synthetic.result[[1]]$V1))),abs(diff(range(synthetic.result[[1]]$V2))))
midx <- sum(range(synthetic.result[[1]]$V1))/2
midy <- sum(range(synthetic.result[[1]]$V2))/2
x <- c(midx-len/2,midx+len/2)*1.1
y <- c(midy-len/2,midy+len/2)*1.1
synthetic.plot[[1]] <- ggplot(data.frame(synthetic.result[[1]]),aes(x=V1,y=V2,size=V3))+
  geom_jitter(color=color.synthetic,width = 2, height=2)+
  ggtitle(plot.title[1])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[1],1))+ylab(paste(plot.title[1],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=20))+theme(legend.position='none')+
  coord_equal(xlim=x,ylim=y)
for (i in 2:9){
  len <- max(abs(diff(range(synthetic.result[[i]]$V1))),abs(diff(range(synthetic.result[[i]]$V2))))
  midx <- sum(range(synthetic.result[[i]]$V1))/2
  midy <- sum(range(synthetic.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  synthetic.plot[[i]] <- ggplot(data.frame(synthetic.result[[i]]),aes(x=V1,y=V2))+
    geom_point(color=color.synthetic[synthetic.label],size=0.8)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}

###################### PBMC3K #######################
pbmc.3k <- Read10X(data.dir= paste0(dir,"/data/pbmc3k/filtered_gene_bc_matrices/hg19/"))
pbmc <- RunSeurat(pbmc.3k,npcs=30)
# Annotate Celltypes
new.cluster.ids <- c("Naive CD4 T","Memory CD4 T","CD14+ Mono","B","FCGR3A+Mono","Effector CD8 T",
                     "NK","Memory CD8 T","Naive CD8 T","DC","IFN-activated CD4T","Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc.label    <- as.data.frame(pbmc@active.ident)
pbmc.label    <- cbind(pbmc.label,as.numeric(pbmc@active.ident))
pbmc.pca      <- pbmc@reductions$pca@cell.embeddings
pbmc.umap     <- pbmc@reductions$umap@cell.embeddings
pbmc.tsne     <- pbmc@reductions$tsne@cell.embeddings
pbmc.supumap <- read.csv(paste0(result.dir,'pbmc3k_supUMAP.csv'),header=F)
pbmc.suppca <- read.csv(paste0(result.dir,'pbmc3k_supPCA.csv'),header=F)
pbmc.mds <- read.csv(paste0(result.dir,'pbmc3k_MDS.csv'),header=F)
pbmc.cpm      <- supCPM(pbmc.pca[,1:30],pbmc.label[,2],dist='geodesic',degree=1,factor=1,niter2=0,init=F)
pbmc.supcpm_euclidean  <- supCPM(pbmc.pca[,1:30],pbmc.label[,2],ratio=0.9,seed=40,init=F,niter2=1000,niter1=400)
pbmc.supcpm_geodesic <- supCPM(pbmc.pca[,1:30],pbmc.label[,2],ratio=0.9,seed=40,dist='geodesic',degree=1,init=F,niter2=1000,niter1=400)

## metrics
for(i in 2:9){
  metric = eval(str2lang(paste0("pbmc.",tolower(plot.title[i]))))%>%VisualMetric(pbmc.pca[,1:30],pbmc.label[,2],4,10)
  for(j in 1:5){
    metric.data.frame[k,1:3] = c('PBMC',plot.title[[i]],names(metric)[[j]])
    metric.data.frame[k,4] = metric[[j]]
    k<-k+1
  }
}


## Plot
pbmc.celltype <- c("Naive CD4 T","Memory CD4 T","CD14+ Mono","B","FCGR3A+Mono","Effector CD8 T","NK",
                   "Memory CD8 T","Naive CD8 T","DC","IFN-activated CD4T","Platelet")
color.pbmc <- rgb(matrix(c(240,128,128, 178,34,34, 0,206,209,
                           50,205,50,   65,105,225,255,215,0,
                           105,105,105, 255,140,0, 240,230,140,
                           147,112,219, 255,69,0,  255,105,180),ncol=3,byrow = T), maxColorValue = 255)
pbmc3k.result <- list()
for(i in 1:9){
  result <- eval(str2lang(paste0("pbmc.",tolower(plot.title[i]))))
  colnames(result) <- paste0('V',1:ncol(result))
  pbmc3k.result[[i]] <- as.data.frame(result)
}
pbmc.plot <- list()
len <- max(abs(diff(range(pbmc3k.result[[1]]$V1))),abs(diff(range(pbmc3k.result[[1]]$V2))))
midx <- sum(range(pbmc3k.result[[1]]$V1))/2
midy <- sum(range(pbmc3k.result[[1]]$V2))/2
x <- c(midx-len/2,midx+len/2)
y <- c(midy-len/2,midy+len/2)
pbmc.plot[[1]] <- ggplot(data.frame(pbmc3k.result[[1]][12:1,]),aes(x=V1,y=V2,size=V3,color=pbmc.celltype[12:1]))+
  geom_point()+
  ggtitle(plot.title[1])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[1],1))+ylab(paste(plot.title[1],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=20))+labs(color ="Celltypes")+
  scale_colour_manual(values=color.pbmc[c(4,3,10,6,5,11,2,8,1,9,7,12)])+
  theme(legend.position = c(0.75, 0.2),legend.key.size = unit(13,'pt'),
        legend.title = element_text(size=8))+
  guides(size=FALSE,color=guide_legend(ncol=2))+
  coord_equal(xlim=x,ylim=y)
for (i in 2:9){
  len <- max(abs(diff(range(pbmc3k.result[[i]]$V1))),abs(diff(range(pbmc3k.result[[i]]$V2))))
  midx <- sum(range(pbmc3k.result[[i]]$V1))/2
  midy <- sum(range(pbmc3k.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  pbmc.plot[[i]] <- ggplot(data.frame(pbmc3k.result[[i]]),aes(x=V1,y=V2))+
    geom_point(color=color.pbmc[pbmc.label[,2]],size=0.5)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}
# add highlight
pbmc.plot[[1]] <- pbmc.plot[[1]]+
  annotate("rect", xmin=3472-150, xmax=3472+150, ymin=153-150, ymax=153+150,
            fill=NA, colour="black")
pbmc.plot[[2]] <- pbmc.plot[[2]]+
  annotate("rect", xmin=5.06-3, xmax=5.06+3, ymin=-8-3, ymax=-8+3,
            fill=NA, colour="black")
pbmc.plot[[3]] <- pbmc.plot[[3]]+
  annotate("rect", xmin=-4.8-3, xmax=-4.8+3, ymin=2.66-3, ymax=2.66+3,
           fill=NA, colour="black")
pbmc.plot[[4]] <- pbmc.plot[[4]]+
  annotate("rect", xmin=10-2, xmax=10+2, ymin=11-2, ymax=11+2,
           fill=NA, colour="black")
pbmc.plot[[5]] <- pbmc.plot[[5]]+
  annotate("rect", xmin=6.5-1, xmax=6.5+1, ymin=-11-1, ymax=-11+1,
           fill=NA, colour="black")
pbmc.plot[[6]] <- pbmc.plot[[6]]+
  annotate("rect", xmin=1.95-1, xmax=1.95+1, ymin=16.5-1, ymax=16.5+1,
           fill=NA, colour="black")
pbmc.plot[[7]] <- pbmc.plot[[7]]+
  annotate("rect", xmin=7.1-2, xmax=7.1+2, ymin=2.56-2, ymax=2.56+2,
           fill=NA, colour="black")
pbmc.plot[[8]] <- pbmc.plot[[8]]+
  annotate("rect", xmin=13.5-2, xmax=13.5+1, ymin=14-1, ymax=14+1,
           fill=NA, colour="black")
pbmc.plot[[9]] <- pbmc.plot[[9]]+
  annotate("rect", xmin=8.5-0.6, xmax=8.5+0.6, ymin=5.6-1, ymax=5.6+0.5,
           fill=NA, colour="black")

###################### RNAmix #######################
rnamix.data <- read.csv(paste0(dir,"/data/RNAmix/RNAmix_celseq2.count.csv"))
rnamix <- RunSeurat(rnamix.data, npcs = 10, norm='log', resolution = 1, nfeature=2000)
rnamix.label <- as.numeric(rnamix@active.ident)
rnamix.pca       <- rnamix@reductions$pca@cell.embeddings
rnamix.umap      <- rnamix@reductions$umap@cell.embeddings
rnamix.tsne      <- rnamix@reductions$tsne@cell.embeddings
rnamix.sumap     <- umap(rnamix.pca[,1:10],y=rnamix.label,metric='cosine')
rnamix.cpm       <- supCPM(rnamix.pca[,1:10],rnamix.label,dist='geodesic',degree=1,factor=1,niter2=0,init=F)
rnamix.supcpm_euclidean   <- supCPM(rnamix.pca[,1:10],rnamix.label,niter2=600,ratio=0.8,init=T)
rnamix.supcpm_geodesic  <- supCPM(rnamix.pca[,1:10],rnamix.label,ratio=0.7,dist='geodesic',degree=1,init=T)
rnamix.supumap <- read.csv(paste0(result.dir,'RNAmix_supUMAP.csv'),header=F)
rnamix.suppca <- read.csv(paste0(result.dir,'RNAmix_supPCA.csv'),header=F)
rnamix.mds <- read.csv(paste0(result.dir,'RNAmix_MDS.csv'),header=F)


## metrics
for(i in 2:9){
  metric = eval(str2lang(paste0("rnamix.",tolower(plot.title[i]))))%>%VisualMetric(rnamix.pca[,1:10],rnamix.label,3,10)
  for(j in 1:5){
    metric.data.frame[k,1:3] = c('RNAmix',plot.title[[i]],names(metric)[[j]])
    metric.data.frame[k,4] = metric[[j]]
    k<-k+1
  }
}


## Plot
color.RNAmix <- rgb(matrix(c(118,42,131,102,189,99,244,109,67,215,48,39,116,173,209,33,102,172,27,120,55),ncol=3,byrow = T),maxColorValue = 255)
RNAmix.result <- list()
for(i in 1:9){
  result <- eval(str2lang(paste0("rnamix.",tolower(plot.title[i]))))
  colnames(result) <- paste0('V',1:ncol(result))
  RNAmix.result[[i]] <- as.data.frame(result)
}

# Plot for RNAmix

RNAmix.plot <- list()
# range for x-axis and y-axis
len <- max(abs(diff(range(RNAmix.result[[1]]$V1))),abs(diff(range(RNAmix.result[[1]]$V2))))
midx <- sum(range(RNAmix.result[[1]]$V1))/2
midy <- sum(range(RNAmix.result[[1]]$V2))/2
x <- c(midx-len/2,midx+len/2)
y <- c(midy-len/2,midy+len/2)
RNAmix.plot[[1]] <- ggplot(data.frame(RNAmix.result[[1]]),aes(x=V1,y=V2,size=V3))+
  geom_point(color=color.RNAmix)+
  ggtitle(plot.title[1])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[1],1))+ylab(paste(plot.title[1],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=20))+theme(legend.position='none')+
  coord_equal(xlim=x,ylim=y)
for (i in 2:9){
  len <- max(abs(diff(range(RNAmix.result[[i]]$V1))),abs(diff(range(RNAmix.result[[i]]$V2))))
  midx <- sum(range(RNAmix.result[[i]]$V1))/2
  midy <- sum(range(RNAmix.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  RNAmix.plot[[i]] <- ggplot(data.frame(RNAmix.result[[i]]),aes(x=V1,y=V2))+
    geom_point(color=color.RNAmix[rnamix.label],size=0.8)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}


##################### Cancer #######################
cancer_data <- load_mixseq_data(paste0(dir,"/SupCPM/data/rna/DMSO_24hr_expt1"),"cancer")
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
# Run results
cancer.label     <- data.frame(as.factor(cancer@active.ident),as.numeric(cancer@active.ident))
cancer.pca       <- cancer@reductions$pca@cell.embeddings
cancer.umap      <- cancer@reductions$umap@cell.embeddings
cancer.tsne      <- cancer@reductions$tsne@cell.embeddings
cancer.cpm       <- supCPM(cancer.pca[,1:20],cancer.label[,2],dist='geodesic',degree=1,k=30,factor=1,niter2=0,init=F)
cancer.supcpm_euclidean   <- supCPM(cancer.pca[,1:20],cancer.label[,2],ratio=0.6,init=F)
cancer.supcpm_geodesic  <- supCPM(cancer.pca[,1:20],cancer.label[,2],ratio=0.6,dist='geodesic',degree=1,init=F,k=30)
cancer.supumap <- read.csv(paste0(result.dir,'cancer_supUMAP.csv'),header=F)
cancer.suppca <- read.csv(paste0(result.dir,'cancer_supPCA.csv'),header=F)
cancer.mds <- read.csv(paste0(result.dir,'cancer_MDS.csv'),header=F)



## metrics
for(i in 2:9){
  metric = eval(str2lang(paste0("cancer.",tolower(plot.title[i]))))%>%VisualMetric(cancer.pca[,1:20],cancer.label[,2],4,10)
  for(j in 1:5){
    metric.data.frame[k,1:3] = c('Cancer',plot.title[[i]],names(metric)[[j]])
    metric.data.frame[k,4] = metric[[j]]
    k<-k+1
  }
}

## Plot
color.cancer <- rgb(matrix(c(97,156,255, 159,71,36,   0,0,0,
                             152,174,12, 181,230,29,  255,201,14,
                             0,186,56,   219,114,251, 252,234,67,
                             21,191,229, 253,98,195),ncol=3,byrow = T),maxColorValue = 255)
cancer.result <- list()
for(i in 1:9){
  result <- eval(str2lang(paste0("cancer.",tolower(plot.title[i]))))
  colnames(result) <- paste0('V',1:ncol(result))
  cancer.result[[i]] <- as.data.frame(result)
}


cancer.celltype <- c("Upper Aerodigestive Tract 1", "Kidney", "Prostate",
                     "Lung 1", "Lung 2", "Nervous System 1",
                     "Lung 3", "Large Intestine 1", "Nervous System 2",
                     "Upper Aerodigestive Tract 2", "Large Intestine 2")

# Plot for cancer
cancer.plot <- list()
len <- max(abs(diff(range(cancer.result[[1]]$V1))),abs(diff(range(cancer.result[[1]]$V2))))
midx <- sum(range(cancer.result[[1]]$V1))/2
midy <- sum(range(cancer.result[[1]]$V2))/2
x <- c(midx-len/2,midx+len/2)
y <- c(midy-len/2,midy+len/2)
cancer.plot[[1]] <- ggplot(data.frame(cancer.result[[1]]),aes(x=V1,y=V2,size=V3,color=cancer.celltype))+geom_point()+
  ggtitle(plot.title[1])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[1],1))+ylab(paste(plot.title[1],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=20))+labs(color ="Celltypes")+
  scale_colour_manual(values=color.cancer[c(2,8,11,4,5,7,6,9,3,1,10)])+
  theme(legend.position = c(0.83, 0.70),legend.key.size = unit(13,'pt'))+guides(size="none")+
  coord_equal(xlim=x,ylim=y)
for (i in 2:9){
  len <- max(abs(diff(range(cancer.result[[i]]$V1))),abs(diff(range(cancer.result[[i]]$V2))))
  midx <- sum(range(cancer.result[[i]]$V1))/2
  midy <- sum(range(cancer.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  cancer.plot[[i]] <- ggplot(data.frame(cancer.result[[i]]),aes(x=V1,y=V2))+
    geom_point(color=color.cancer[cancer.label[,2]],size=0.7)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}

####################### COVID #######################
meta.t <- read.table(paste0(dir,"/COVID/NKT.cell.annotation.meta.txt"),header=T,sep=';')
nCoV <- readRDS(paste0(dir,"/COVID/nCoV.rds"))
meta.t <- meta.t[-which(meta.t$celltype=='Doublets'|meta.t$celltype=='Uncertain'),]
label.t <- meta.t$celltype
NKT.label <- as.numeric(as.factor(label.t))
id <- which(colnames(nCoV)%in%meta.t$ID)
CovidT <- nCoV[,id]
CovidT <- RunPCA(CovidT,verbose = F)
CovidT <- RunTSNE(CovidT,dims=1:30)
CovidT <- RunUMAP(CovidT,dims=1:30)
CovidT.pca  <- CovidT@reductions$pca@cell.embeddings
CovidT.tsne <- CovidT@reductions$tsne@cell.embeddings
CovidT.umap <- CovidT@reductions$umap@cell.embeddings
CovidT.euclidean <- supCPM(CovidT.pca[,1:30],NKT.label,dist='euclidean',niter1=300,niter2=900,seed=40,ratio=0.85,factor=1.3,degree=3,init=F)
CovidT.geodesic <- supCPM(CovidT.pca[,1:30],NKT.label,dist='geodesic',niter1=300,niter2=900,seed=40,ratio=0.85,degree=2,init=F)
CovidT.supumap <- read.csv(paste0(result.dir,'CovidT_supUMAP.csv'),header=F)
CovidT.suppca <- read.csv(paste0(result.dir,'CovidT_supPCA.csv'),header=F)
CovidT.mds <- read.csv(paste0(result.dir,'CovidT_MDS.csv'),header=F)

for(i in 2:9){
  metric = eval(str2lang(paste0("CovidT.",tolower(plot.title[i]))))%>%VisualMetric(nT.pca[,1:30],NKT.label,3,10)
  for(j in 1:5){
    metric.data.frame[k,1:3] = c('CovidT',plot.title[[i]],names(metric)[[j]])
    metric.data.frame[k,4] = metric[[j]]
    k<-k+1
  }
}


CovidT.celltype <- c("Tfh","Tim3-CD8+IFN-gamma(low) T-cells","Th17","gamma-delta T-cells",
                  "Tim3+CD8+IFN-gamma(high) T-cells","T reg")
CovidT.color <- rgb(matrix(c(223,111,103,230,143,22,50,147,255,
                          70,200,98,150,161,48,0,171,190),ncol=3,byrow = T),maxColorValue = 255)
CovidT.result <- list()
for(i in 1:9){
  result <- eval(str2lang(paste0("CovidT.",tolower(plot.title[i]))))
  colnames(result) <- paste0('V',1:ncol(result))
  CovidT.result[[i]] <- as.data.frame(result)
}

CovidT.plot <- list()
len <- max(abs(diff(range(CovidT.result[[1]]$V1))),abs(diff(range(CovidT.result[[1]]$V2))))
midx <- sum(range(CovidT.result[[1]]$V1))/2
midy <- sum(range(CovidT.result[[1]]$V2))/2
x <- c(midx-len/2,midx+len/2)
y <- c(midy-len/2,midy+len/2)
CovidT.plot[[1]] <- ggplot(data.frame(CovidT.result[[1]]))+
  geom_point(aes(x=V1,y=V2,size=V3,color=CovidT.celltype))+
  ggtitle(plot.title[1])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[1],1))+ylab(paste(plot.title[1],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=20))+labs(color ="Celltypes")+
  scale_color_manual(values=CovidT.color)+
  theme(legend.position = c(0.51, 0.15),legend.key.size = unit(15,'pt'),
        legend.title = element_text(size=11))+
  guides(size=FALSE,color=guide_legend(ncol=2))+
  coord_equal(xlim=x,ylim=y)
CovidT.plot[[1]]
for (i in 2:9){
  len <- max(abs(diff(range(CovidT.result[[i]]$V1))),abs(diff(range(CovidT.result[[i]]$V2))))
  midx <- sum(range(CovidT.result[[i]]$V1))/2
  midy <- sum(range(CovidT.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  CovidT.plot[[i]] <- ggplot(data.frame(CovidT.result[[i]]),aes(x=V1,y=V2))+
    geom_point(color=CovidT.color[CovidT.label],size=0.5)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)+theme(legend.position="none")
}

##
metric.data.frame <- read.csv("C:/SupCPM/Metric.csv")
metric.data.frame <- metric.data.frame[,-1]
metric.data.frame[metric.data.frame$dataset=='NKT',1]='CovidT'
metric.data.frame.part <- metric.data.frame %>% filter(method!='supPCA'&method!='PCA')
metric.data.frame.part$dataset <- factor(metric.data.frame.part$dataset,levels=c("Synthetic","RNAmix","Cancer","PBMC","CovidT"))
levels(metric.data.frame.part$dataset) <- c("Synthetic","RNAmix","Cancer","PBMC","CovidT")
metric.data.frame.part[metric.data.frame.part$metric=='CS','Scores'] =  -0.5*log(metric.data.frame.part[metric.data.frame.part$metric=='CS','Scores'])
metric.data.frame.part[metric.data.frame.part$metric=='CS','metric'] = '-0.5log(CSM)'
metric.plot <- ggplot(metric.data.frame.part,aes(x=dataset,y=Scores,fill=method))+
  geom_bar(position = 'dodge',stat='identity')+
  facet_wrap(~metric)+
  theme(strip.text.x = element_text(size=16,face = "bold.italic"),axis.title.x = element_blank(),axis.text.x = element_text(size=15),axis.title.y = element_text(size=13))+
  scale_fill_brewer(palette="Set2")+theme(legend.key.size = unit(1.5, 'cm'))+
  theme(legend.title = element_text(size=20))+
  theme(legend.text = element_text(size=18))


## save the plots
pdf(file = paste0(fig.dir,"Cancer Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = cancer.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file = paste0(fig.dir,"PBMC Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = pbmc.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file = paste0(fig.dir,"RNAmix Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = RNAmix.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file = paste0(fig.dir,"Synthetic Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = synthetic.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file = paste0(fig.dir,"CovidT Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = CovidT.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file = paste0(fig.dir,"Metric.pdf"),width=25,height=15)
print(metric.plot)
dev.off()


######################  EB  ########################
library(reticulate)
source_python("EB.py")
EBdata <- EB_preprocess(paste0(dir,"/COVID"))
count <- EBdata[[1]]
label <- EBdata[[2]]
rm(EBdata)
row.names(count) <- 1:nrow(count)
colnames(count) <- 1:ncol(count)
set.seed(123)
id <- sample(1:nrow(count),3000)
EBpart <- count[id,]
labelEB <- label[id]
EB.label <- as.numeric(as.factor(labelEB))
#labelEB.refine <- label.refine[id]
EB <- RunSeurat(t(EBpart),npcs=20,nfeature = 3000)
EB.pca      <- EB@reductions$pca@cell.embeddings
EB.umap     <- EB@reductions$umap@cell.embeddings
EB.tsne     <- EB@reductions$tsne@cell.embeddings
EB.scale    <- as.matrix(EB@assays$SCT@scale.data)
EB.mds <- read.csv(paste0(result.dir,'EB_MDS.csv'),header=F)
EB.suppca <- read.csv(paste0(result.dir,'EB_supPCA.csv'),header=F)
EB.supumap <- read.csv(paste0(result.dir,'EB_supUMAP.csv'),header=F)
EB.supcpm_euclidean <- supCPM(EB.pca[,1:20],EB.label,dist='euclidean',niter1=300,niter2=700,seed=40,ratio=0.7,factor=1.3,init=F)
EB.supcpm_geodesic <- supCPM(EB.pca[,1:20],EB.label,dist='geodesic',k=10, niter1=300,niter2=700,degree=1,seed=40,ratio=0.7,init=F)
EB.cpm <- supCPM(EB.pca[,1:20],EB.label,dist='geodesic',k=10, niter1=300,niter2=0,degree=1,seed=40,ratio=0,init=F,factor=1)

EB.celltype <- factor(c("Day 0-3","Day 6-9","Day 12-15","Day 18-21","Day 24-27"))
levels(EB.celltype) <- c("Day 0-3","Day 6-9","Day 12-15","Day 18-21","Day 24-27")
EB.result <- list()
for(i in 1:9){
  result <- eval(str2lang(paste0("EB.",tolower(plot.title[i]))))
  colnames(result) <- paste0('V',1:ncol(result))
  EB.result[[i]] <- as.data.frame(result)
}


EB.plot <- list()
len <- max(abs(diff(range(EB.result[[1]]$V1))),abs(diff(range(EB.result[[1]]$V2))))
midx <- sum(range(EB.result[[1]]$V1))/2
midy <- sum(range(EB.result[[1]]$V2))/2
x <- c(midx-len/2,midx+len/2)
y <- c(midy-len/2,midy+len/2)
EB.plot[[1]] <- ggplot(data.frame(EB.result[[1]][c(1,5,2,3,4),]),aes(x=V1,y=V2,size=V3,color=EB.celltype))+
  geom_point()+
  ggtitle(plot.title[1])+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[1],1))+ylab(paste(plot.title[1],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=20))+labs(color ="Time")+
  scale_color_viridis_d()+
  theme(legend.position = c(0.3, 0.2),legend.key.size = unit(15,'pt'),
        legend.title = element_text(size=11))+
  guides(size=FALSE,color=guide_legend(ncol=2))+
  coord_equal(xlim=x,ylim=y)
for (i in 2:9){
  len <- max(abs(diff(range(EB.result[[i]]$V1))),abs(diff(range(EB.result[[i]]$V2))))
  midx <- sum(range(EB.result[[i]]$V1))/2
  midy <- sum(range(EB.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  EB.plot[[i]] <- ggplot(data.frame(EB.result[[i]]),aes(x=V1,y=V2))+
    geom_point(aes(color=as.factor(EB.label)),size=0.5)+
    ggtitle(plot.title[i])+
    scale_color_viridis_d()+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)+theme(legend.position="none")
}

pdf(file = paste0(fig.dir,"EB Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = EB.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()













