## Plot the result
library(ggplot2)
library(cowplot)
output.dir = '.../result/output/'
figure.dir = '.../result/figures/'
# color
color.pbmc <- rgb(matrix(c(240,128,128, 178,34,34, 0,206,209,
                           50,205,50,   65,105,225,255,215,0, 
                           105,105,105, 255,140,0, 240,230,140,
                           147,112,219, 255,69,0,  255,105,180),ncol=3,byrow = T), maxColorValue = 255)
color.RNAmix <- rgb(matrix(c(118,42,131,102,189,99,244,109,67,215,48,39,116,173,209,33,102,172,27,120,55),ncol=3,byrow = T),maxColorValue = 255)
color.synthetic <- rgb(matrix(c(0,136,64,99,189,95,167,218,86,255,74,70,
                             163,0,89,255,46,127,27,229,255,0,110,165,
                             0,0,165),ncol=3,byrow = T),maxColorValue = 255)
color.cancer <- rgb(matrix(c(97,156,255, 159,71,36,   0,0,0,
                             152,174,12, 181,230,29,  255,201,14,
                             0,186,56,   219,114,251, 252,234,67,
                             21,191,229, 253,98,195),ncol=3,byrow = T),maxColorValue = 255)

# self-defined functions to load data
load.result <- function(name){
  result.list <- list()
  label <- read.csv(paste0(output.dir,name,'_label.csv'), header = T)
  result.list[[1]] <- read.csv(paste0(output.dir,name,'_MDS.csv'), header = F)
  result.list[[2]] <- read.csv(paste0(output.dir,name,'_pca.csv'), header = T)
  result.list[[2]] <- result.list[[2]][,2:3]
  colnames(result.list[[2]]) <- c('V1','V2')
  result.list[[3]] <- read.csv(paste0(output.dir,name,'_supPCA.csv'), header = F)
  result.list[[4]] <- read.csv(paste0(output.dir,name,'_tsne.csv'), header = T)
  result.list[[4]] <- result.list[[4]][,-1]
  colnames(result.list[[4]]) <- c('V1','V2')
  result.list[[5]] <- read.csv(paste0(output.dir,name,'_umap.csv'), header = T)
  result.list[[5]] <- result.list[[5]][,-1]
  colnames(result.list[[5]]) <- c('V1','V2')
  result.list[[6]] <- read.csv(paste0(output.dir,name,'_supUMAP.csv'), header = F)
  result.list[[7]] <- read.csv(paste0(output.dir,name,'_CPM.csv'), header = F)
  result.list[[8]] <- read.csv(paste0(output.dir,name,'_Eucli.csv'), header = F)
  result.list[[9]] <- read.csv(paste0(output.dir,name,'_Geo.csv'), header = F)
  return(list(result.list,label))
}

# load data
RNAmix.result <- load.result('RNAmix')
RNAmix.label <- RNAmix.result[[2]]
RNAmix.result <- RNAmix.result[[1]]

pbmc3k.result <- load.result('pbmc3k')
pbmc3k.label <- pbmc3k.result[[2]]
pbmc3k.result <- pbmc3k.result[[1]]

cancer.result <- load.result('cancer')
cancer.label <- cancer.result[[2]]
cancer.result <- cancer.result[[1]]

synthetic.result <- load.result('synthetic')
synthetic.label <- synthetic.result[[2]]
synthetic.result <- synthetic.result[[1]]

# title & legend
plot.title <- c('MDS','PCA','supPCA','t-SNE','UMAP','supUMAP','CPM','supCPM-Euclidean','supCPM-Geodesic')
pbmc.celltype <- c("Naive CD4 T","CD14+ Mono","Memory CD4 T","B","Effector CD8 T","NK","FCGR3A+Mono",
                   "Memory CD8 T","Naive CD8 T","DC","IFN-activated CD4T","Platelet"
                   )
pbmc.celltype <- c("Naive CD4 T","Memory CD4 T","CD14+ Mono","B","FCGR3A+Mono","Effector CD8 T","NK",
                   "Memory CD8 T","Naive CD8 T","DC","IFN-activated CD4T","Platelet")
cancer.celltype <- c("Upper Aerodigestive Tract 1", "Kidney", "Prostate",
                   "Lung 1", "Lung 2", "Nervous System 1",
                   "Lung 3", "Large Intestine 1", "Nervous System 2",
                   "Upper Aerodigestive Tract 2", "Large Intestine 2")

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
    geom_point(color=color.RNAmix[RNAmix.label[,2]+1],size=0.8)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}

# Plot for pbmc3k
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
    geom_point(color=color.pbmc[pbmc3k.label[,3]],size=0.5)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}
# add highlight
pbmc.plot[[8]] <- pbmc.plot[[8]]+ 
  annotate("rect", xmin=16.5, xmax=18.5, ymin=-10.2, ymax=-8.8, 
           fill=NA, colour="black")
pbmc.plot[[9]] <- pbmc.plot[[9]]+ 
  annotate("rect", xmin=1, xmax=2.5, ymin=-10, ymax=-8.5, 
           fill=NA, colour="black")
pbmc.plot[[4]] <- pbmc.plot[[4]]+
  annotate("rect", xmin=8, xmax=11.5, ymin=10, ymax=13, 
           fill=NA, colour="black")
pbmc.plot[[5]] <- pbmc.plot[[5]]+
  annotate("rect", xmin=5.5, xmax=7.5, ymin=-11.5, ymax=-9.5, 
           fill=NA, colour="black")
pbmc.plot[[6]] <- pbmc.plot[[6]]+
  annotate("rect", xmin=0.5, xmax=3.5, ymin=15, ymax=18, 
           fill=NA, colour="black")

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
  theme(legend.position = c(0.83, 0.70),legend.key.size = unit(13,'pt'))+guides(size=FALSE)+
  coord_equal(xlim=x,ylim=y)
for (i in 2:9){
  len <- max(abs(diff(range(cancer.result[[i]]$V1))),abs(diff(range(cancer.result[[i]]$V2))))
  midx <- sum(range(cancer.result[[i]]$V1))/2
  midy <- sum(range(cancer.result[[i]]$V2))/2
  x <- c(midx-len/2,midx+len/2)
  y <- c(midy-len/2,midy+len/2)
  cancer.plot[[i]] <- ggplot(data.frame(cancer.result[[i]]),aes(x=V1,y=V2))+
    geom_point(color=color.cancer[cancer.label[,3]],size=0.7)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
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
    geom_point(color=color.synthetic[synthetic.label[,2]],size=0.8)+
    ggtitle(plot.title[i])+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
    xlab(paste(plot.title[i],1))+ylab(paste(plot.title[i],2))+
    scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
    theme(plot.title = element_text(hjust = 0.5,size=20))+
    coord_equal(xlim=x,ylim=y)
}


# save
pdf(file = paste0(figure.dir,"Synthetic Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = synthetic.plot,ncol=2,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95,axis='l')
print(a)
dev.off()

pdf(file = paste0(figure.dir,"PBMC3k Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = pbmc.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file =  paste0(figure.dir,"RNAmix Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = RNAmix.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()

pdf(file = paste0(figure.dir,"Cancer Data.pdf"),width=15,height=15)
a=plot_grid(plotlist = cancer.plot,ncol=3,labels = letters[1:9],
            label_size = 20,label_fontface ='bold',scale=0.95)
print(a)
dev.off()



## Plot for introduction
len <- max(abs(diff(range(RNAmix.result[[5]]$V1))),abs(diff(range(RNAmix.result[[5]]$V2))))
midx <- sum(range(RNAmix.result[[5]]$V1))/2
midy <- sum(range(RNAmix.result[[5]]$V2))/2
x <- c(midx-len/2,midx+len/2)
y <- c(midy-len/2,midy+len/2)
umap <- ggplot(data.frame(RNAmix.result[[5]]),aes(x=V1,y=V2))+
  geom_point(color=color.RNAmix[RNAmix.label[,2]+1],size=2.5)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_line(colour='black'))+
  xlab(paste(plot.title[5],1))+ylab(paste(plot.title[5],2))+
  scale_x_continuous(labels = NULL)+scale_y_continuous(labels = NULL)+
  theme(plot.title = element_text(hjust = 0.5,size=12),axis.title = element_text(size=20))+
  coord_equal(xlim=x,ylim=y)+annotate("rect", xmin=-0.5, xmax=-0, ymin=0.3, ymax=0.7, 
                                      fill=NA, colour="black")+
  annotate("rect", xmin=-0.25, xmax=0.25, ymin=-1.2, ymax=-0.83, 
            fill=NA, colour="black")
pdf(file = paste0(figure.dir,"introduction.pdf"),width=10,height=10)
print(umap)
dev.off()


## Metric Plot
method <- c('supCPM-Geo','supCPM-Eu','CPM','t-SNE','supUMAP','PCA','supPCA','UMAP')
metric <- c('CorV','KNC','CPD','logFisher','KNN')
dataset <- c('PBMC','RNAmix','Cancer Cell Lines','Synthetic')
metrics.rnamix <- read.table(paste0(output.dir,'RNAmix_metric.csv'),sep=',')
metrics.cancer <- read.table(paste0(output.dir,'cancer_metric.csv'),sep=',')
metrics.pbmc <- read.table(paste0(output.dir,'pbmc3k_metric.csv'),sep=',')
metrics.synthetic <- read.table(paste0(output.dir,'synthetic_metric.csv'),sep=',')

metric.data.frame <- data.frame(matrix(rep(0,160*4),ncol=4))
colnames(metric.data.frame) <- c('dataset','method','metric','Scores')
metric.data.frame[,1] <- rep(dataset,each=40)
metric.data.frame[,2] <- rep(method,each=5)
metric.data.frame[,3] <- rep(metric,32)
metric.data.frame[1:40,4] <- as.vector(as.matrix(metrics.pbmc))
metric.data.frame[41:80,4] <- as.vector(as.matrix(metrics.rnamix))
metric.data.frame[81:120,4] <- as.vector(as.matrix(metrics.cancer))
metric.data.frame[121:160,4] <- as.vector(as.matrix(metrics.synthetic))
metric.data.frame[metric.data.frame$metric=='logFisher','Scores'] <--0.1*log(metric.data.frame[metric.data.frame$metric=='logFisher','Scores'])
metric.data.frame$method <- factor(metric.data.frame$method, levels=c('CPM','supCPM-Eu','supCPM-Geo','PCA','supPCA','t-SNE','UMAP','supUMAP'), ordered=TRUE)
metric.data.frame$dataset <- factor(metric.data.frame$dataset, levels=c('Synthetic','RNAmix','Cancer Cell Lines','PBMC'), ordered=TRUE)


# need jittering
set.seed(321)
jitter <- c(6,74,85,95,117,120)
not.jitter <- (1:160)[-jitter]
metric.plot <- ggplot(data=metric.data.frame[not.jitter,])+
  geom_point(aes(x=method,y=Scores,color=metric),size=2.2,pch=16,alpha=0.7)+
  geom_jitter(aes(x=method,y=Scores,color=metric),data=metric.data.frame[jitter[1],],size=2.2,pch=16,alpha=0.7,height = 0,width=0.25)+
  geom_jitter(aes(x=method,y=Scores,color=metric),data=metric.data.frame[jitter[2],],size=2.2,pch=16,alpha=0.7,height = 0,width=0.25)+
  geom_jitter(aes(x=method,y=Scores,color=metric),data=metric.data.frame[jitter[3],],size=2.2,pch=16,alpha=0.7,height = 0,width=0.25)+
  geom_jitter(aes(x=method,y=Scores,color=metric),data=metric.data.frame[jitter[4],],size=2.2,pch=16,alpha=0.7,height = 0,width=0.25)+
  geom_jitter(aes(x=method,y=Scores,color=metric),data=metric.data.frame[jitter[5],],size=2.2,pch=16,alpha=0.7,height = 0,width=0.55)+
  geom_jitter(aes(x=method,y=Scores,color=metric),data=metric.data.frame[jitter[6],],size=2.2,pch=16,alpha=0.7,height = 0,width=0.25)+
  facet_wrap(~dataset,ncol=4)+theme(axis.text.x = element_text(angle = 40, hjust = 0.7, vjust = 0.8,size=10))+
  xlab('')+scale_color_manual(values=c('#4C9E00','#E67350','#7E92E6','#CFC03E','#AA68AC'))+
  theme(plot.title = element_text(hjust = 0.5,size=20))

metric.plot
pdf(file = paste0(figure.dir,"Metrics.pdf"),width=12,height=5)
print(metric.plot)
dev.off()
















