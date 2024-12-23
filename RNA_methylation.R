gene<-read.table("gene.txt",head=F,sep='\t',check.names = F)
dat<-read.table("uniq.symbol.txt",head=T,sep='\t',check.names = F,row.names = 1)
expr<-dat[intersect(gene$V1,row.names(dat)),]
data<-as.matrix(expr)
maxK=9
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(data,
                               maxK=maxK,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title="Cluster",
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123456,
                               plot="pdf")
Kvec = 2:9
x1 = 0.1; x2 = 0.9 
PAC = rep(NA,length(Kvec)) 
names(PAC) = paste("K=",Kvec,sep="") 
for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
}
optK = Kvec[which.min(PAC)]
optK
PAC <- as.data.frame(PAC)
PAC$K <- 2:9
library(ggplot2)
ggplot(PAC,aes(factor(K),PAC,group=1))+
    geom_line(color="#00A087FF",size=2)+
    theme_bw(base_rect_size = 1.5)+
    geom_point(size=6,shape=21,color="#E64B35FF",fill="#E64B35FF")+
    ggtitle('Proportion of ambiguous clustering')+
    xlab('Cluster number K')+ylab(NULL)+
    theme(axis.text = element_text(size=16,color="black",face="bold"),
          plot.title = element_text(size=16,color="black",face="bold",hjust = 0.5),
          axis.title = element_text(size=16,color="black",face="bold"))
clusterNum=2      
cluster=results[[clusterNum]][["consensusClass"]]
sub <- data.frame(Sample=names(cluster),Cluster=cluster)
sub$Cluster <- paste0('C',sub$Cluster)
table(sub$Cluster)
head(sub)
my <- results[[2]][["ml"]]
library(pheatmap)
rownames(my) <- sub$Sample
colnames(my) <- sub$Sample
library(ggsci)
pheatmap(1-my,show_colnames = F,show_rownames = F,
         treeheight_row = 20,treeheight_col = 20,
         clustering_method = 'complete',
         color = colorRampPalette(c("navy", "white"))(50),
         annotation_names_row = F,annotation_names_col = F,
         annotation_row = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_col = data.frame(Cluster=sub$Cluster,row.names = sub$Sample),
         annotation_colors = list(Cluster=c('C2'="#E64B35FF",'C1'="#4DBBD5FF")),legend_breaks = c(0,0.2,0.4,0.6,0.8,1),legend_labels = c(1,0.8,0.6,0.4,0.2,0))
#Survival analysis
custom_theme <- function() {
    theme_survminer() %+replace%
        theme(panel.grid = element_blank(),
              axis.text =   element_text(color='black',size=15, face = "bold"),
              axis.title = element_text(color='black',size=20, face = "bold"),
              panel.background = element_rect(fill = "#EDF8E9")
        )
}
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
    original.p <- p
    if(is.ggplot(original.p)) list.plots <- list(original.p)
    else if(is.list(original.p)) list.plots <- original.p
    else stop("Can't handle an object of class ", class (original.p))
    .set_font <- function(font){
        font <- ggpubr:::.parse_font(font)
        ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
    }
    for(i in 1:length(list.plots)){
        p <- list.plots[[i]]
        if(is.ggplot(p)){
            if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
            if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
            if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
            if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
            if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
            if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
            if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
            list.plots[[i]] <- p
        }
    }
    if(is.ggplot(original.p)) list.plots[[1]]
    else list.plots
}
library(survival)
library("survminer")
rt=read.table("OS.txt",header=T,sep="\t");
diff=survdiff(Surv(OS.time,OS) ~Group,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(OS.time,OS) ~ Group, data = rt)
p<-ggsurvplot(fit, 
              data=rt,
              conf.int=T,conf.int.style='step', size=1.5,
              pval=paste0 ("P < ",2.2e-16),
              pval.size=12,
              legend.title="Group",
              legend.labs=levels(factor(rt[,"Group"])),
              legend = c(0.9, 0.9),
              font.legend=10,
              xlab="Time(years)",
              break.time.by = 2,
              palette = c("jco"),
              surv.median.line = "hv",
              risk.table=T,ggtheme = custom_theme(),
              cumevents=F,
              risk.table.height=.30)
p$table <- customize_labels(
    p$table,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "purple"), 
    font.caption  = c(14, "plain", "orange"),        
    font.x        = c(20, "bold", "black"),          
    font.y        = c(20, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),font.ytickslab = c(20, "bold", "black")
)
p$plot<-customize_labels(
    p$plot,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "purple"), 
    font.caption  = c(14, "plain", "orange"),        
    font.x        = c(25, "bold", "black"),          
    font.y        = c(25, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),
    font.ytickslab = c(20, "bold", "black")
)
library(ComplexHeatmap)
ee<-read.table("Cluster_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
my<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
ee<-ee[,intersect(row.names(my),colnames(ee))]
identical(row.names(my),colnames(ee))
ee <- t(scale(t(ee)))
ee[ee > 2] <- 2 
ee[ee < -2] <- -2 
my$Status <- factor(my$Status)
my$Cluster <- factor(my$Cluster)
my$Stage <- factor(my$Stage)
my$Gender <- factor(my$Gender,levels = c('Female','Male'))
my$Age <- factor(my$Age,levels = c('<=65','>65'))
my$T=factor(my$T)
my$N=factor(my$N)
my$M=factor(my$M,levels = c('M0','M1','NA'))
Cluster <- c("#00a087","#e64b35","#4dbbd5")
names(Cluster) <- levels(my$Cluster)
Age <- c(pal_nejm(alpha = 0.9)(8)[3],'#CF4E27')
names(Age) <- levels(my$Age)
Gender <- c('#E0864A','rosybrown')
names(Gender) <- levels(my$Gender)
Stage <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(Stage) <- levels(my$Stage)
T <- c('cornsilk','paleturquoise','goldenrod','firebrick')
names(T) <- levels(my$T)
N <- c('paleturquoise','goldenrod','firebrick')
names(N) <- levels(my$N)
M<- c('paleturquoise','firebrick',"white")
names(M) <- levels(my$M)
Status <- c('lavenderblush','slategray')
names(Status) <- levels(my$Status)
col<-list(Cluster,Age,Gender,Stage,T,N,M,Status)


Top = HeatmapAnnotation(Cluster=my$Cluster,
                        Age=my$Age,
                        Gender=my$Gender,
                        Stage= my$Stage,T=my$T,N=my$N,N=my$N,
                        Status = my$Status,
                        annotation_legend_param=list(labels_gp = gpar(fontsize = 10),border = T,title_gp = gpar(fontsize = 10,fontface = "bold"),ncol=1),
                        col=list(Cluster = Cluster,
                                 Age = Age,
                                 Gender = Gender,
                                 Stage= Stage,T=T,N=N,M=M,
                                 Status = Status
                        ),
                        show_annotation_name = TRUE,
                        annotation_name_side="left",
                        annotation_name_gp = gpar(fontsize = 10))


Heatmap(ee,name='Z-score',
        cluster_rows = TRUE,top_annotation = Top,
        col=colorRamp2(c(-2,0,2),c('#1d6fae','white','#ef7700')),
        color_space = "RGB",
        cluster_columns = FALSE,border = T,
        row_order=NULL,
        row_names_side = 'left',
        column_order=NULL,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 9),
        column_split = c(rep(1,356),rep(2,192),rep(3,52)),
        gap = unit(1, "mm"),
        column_title = NULL,
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend = TRUE,
        heatmap_legend_param=list(labels_gp = gpar(fontsize = 13), border = T,
                                  title_gp = gpar(fontsize = 13, fontface = "bold")),
        column_gap = unit(2,'mm'))

#PCA analysis
library(ggplot2)
library(ggh4x)
data<-read.table("Cluster_expr.txt",head=T,sep='\t',check.names = F,row.names = 1)
pca1 <- prcomp(t(data),center = TRUE,scale. = TRUE)
write.table(pca1$x[,c(1,2)],"PCA.txt",quote=F,sep='\t')
dt<-read.table("PCA.txt",head=T,sep='\t',check.names = F,row.names = 1)
mytheme <- theme(panel.grid = element_blank(),
                 panel.background = element_blank(),
                 legend.key = element_blank(),
                 legend.position = "top",
                 axis.line = element_line(colour = "grey30"),axis.text=element_text(color="black",face="bold",size=16),
                 axis.ticks.length = unit(1.8, "mm"),axis.title=element_text(color="black",face="bold",size=16),legend.text = element_text(color="black",face="bold",size=18),
                 ggh4x.axis.ticks.length.minor = rel(0.6))
ggplot(dt,aes(x=PC1,y=PC2,fill=Cluster))+
    stat_centroid(aes(xend = PC1, yend = PC2, colour = Cluster),
                  geom = "segment", crop_other = F,
                  alpha=0.3,size = 1,show.legend = F)+
    geom_point(size=3,alpha=0.7,
               color="white",shape = 21,show.legend = T)+
    scale_color_manual(name="",
                       values = c("#4DBBD5FF","#E64B35FF"))+
    scale_fill_manual(name="",
                      values = c("#4DBBD5FF","#E64B35FF"))+
    scale_x_continuous(expand=expansion(add = c(0.7,0.7)))+
    scale_y_continuous(expand=expansion(add = c(0.5,0.5)))+
    guides(x = "axis_truncated",y = "axis_truncated")+mytheme
	
#WGCNA
library(limma)
library(gplots)
library(WGCNA)
data<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
datTraits<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
selectGenes=names(tail(sort(apply(data,1,sd)), n=round(nrow(data)*0.25)))
data=data[selectGenes,]
datExpr0=t(data)
datExpr0<-datExpr0[intersect(row.names(datTraits),row.names(datExpr0)),]
identical(row.names(datTraits),row.names(datExpr0))
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "01.sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

abline(h = 20000, col="red")
dev.off()
sampleTree2 = hclust(dist(datExpr0), method="average")
traitColors = numbers2colors(datTraits, signed = FALSE)
pdf(file="02.sample_heatmap.pdf", width=12, height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

enableWGCNAThreads()   
powers = c(1:20)       
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5,RsquaredCut = 0.90)
pdf(file="03.scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="#00A087FF");
abline(h=0.90,col="red") 

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="#4DBBD5FF")
dev.off()

sft 
softPower =sft$powerEstimate    
adjacency = adjacency(datExpr0, power = softPower)
softPower

TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="04.gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

minModuleSize = 100     
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="05.Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="06.Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()

moduleColors=dynamicColors
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
select = sample(nGenes, size=1000)     
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method="average")
selectColors = moduleColors[select]
#sizeGrWindow(9,9)
plotDiss=selectTOM^softPower
diag(plotDiss)=NA
myheatcol = colorpanel(250, "red", "orange", "lemonchiffon")    
pdf(file="07.TOMplot.pdf", width=7, height=7)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)
dev.off()


moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="08.Module_trait.pdf", width=6.5, height=5.5)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3.5, 8, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=names(datTraits)
geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

trait="Cluster"
traitColumn=match(trait,traitNames)  
for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
        outPdf=paste("09.", trait, "_", module,".pdf",sep="")
        pdf(file=outPdf,width=7,height=7)
        par(mfrow = c(1,1))
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                           abs(geneTraitSignificance[moduleGenes, traitColumn]),
                           xlab = paste("Module Membership in", module, "module"),
                           ylab = paste("Gene significance for ",trait),
                           main = paste("Module membership vs. gene significance\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
        dev.off()
    }
}

probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                           GSPvalue[, Tra])
    names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                         names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                           MMPvalue[, mod])
    names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                         names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)

gene<-read.table("Gene.txt",head=F,sep='\t',check.names = F)
data<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
time<-read.table("time.txt",head=T,sep='\t',check.names = F,row.names = 1)
expr<-data[intersect(gene$V1,row.names(data)),]
expr<-expr[,intersect(row.names(time),colnames(expr))]
rt<-cbind(time,t(expr))
library(survival)
pFilter=0.05                                                    

outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(coxP<pFilter){
        sigGenes=c(sigGenes,i)
        outTab=rbind(outTab,
                     cbind(id=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
        )
    }
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

work.path <- "E:/RNA methylation"; setwd(work.path) 


code.path <- file.path(work.path, "Codes") 
data.path <- file.path(work.path, "InputData") 
res.path <- file.path(work.path, "Results") 
fig.path <- file.path(work.path, "Figures") 

if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# BiocManager::install("mixOmics")
# BiocManager::install("survcomp")
# devtools::install_github("binderh/CoxBoost")
# install.packages("randomForestSRC")
# install.packages("snowfall")


library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)

source(file.path(code.path, "ML.R"))

FinalModel <- c("panML", "multiCox")[2]

## Training Cohort ---------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) #

# 按队列对数据分别进行标准化（根据情况调整centerFlags和scaleFlags）
## data: 需要表达谱数据（行为样本，列为基因） 
## cohort：样本所属队列，为向量，不输入值时默认全表达矩阵来自同一队列
## centerFlag/scaleFlags：是否将基因均值/标准差标准化为1；
##        默认参数为NULL，表示不进行标准化；
##        为T/F时，表示对所有队列都进行/不进行标准化
##        输入由T/F组成的向量时，按顺序对队列进行处理，向量长度应与队列数一样
##        如centerFlags = c(F, F, F, T, T)，表示对第4、5个队列进行标准化，此时flag顺序应当与队列顺序一致
##        如centerFlags = c("A" = F, "C" = T, "B" = F)，表示对队列C进行标准化，此时不要求flag顺序与data一致
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)
# summary(apply(Train_set, 2, var))
# summary(apply(Test_set, 2, var))
# lapply(split(as.data.frame(Test_set), Test_surv$Cohort), function(x) summary(apply(x, 2, var))) # 测试scale结果

# Model training and validation -------------------------------------------

## method list --------------------------------------------------------
# 此处记录需要运行的模型，格式为：算法1名称[算法参数]+算法2名称[算法参数]
# 目前仅有StepCox和RunEnet支持输入算法参数
methods <- read.xlsx(file.path(code.path, "41467_2022_28421_MOESM4_ESM.xlsx"), startRow = 2)$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------
min.selected.var <- 5 # 筛选变量数目的最小阈值
timeVar = "OS.time"; statusVar = "OS" # 定义需要考虑的结局事件，必须出现在Train_surv以及Test_surv中

## Pre-training 
Variable = colnames(Train_expr)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))
preTrain.method

set.seed(seed = 123) # 设置建模种子，使得结果可重复
preTrain.var <- list()
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 机器学习方法
                                 Train_expr = Train_set, # 训练集有潜在预测价值的变量
                                 Train_surv = Train_surv, # 训练集生存数据
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的生存变量，必须出现在Train_surv中
}
preTrain.var[["simple"]] <- colnames(Train_expr)

model <- list() # 初始化模型结果列表
set.seed(seed = 123) # 设置建模种子，使得结果可重复
for (method in methods){ # 循环每一种方法组合
  # method <- "CoxBoost+plsRcox" # [举例]若遇到报错，请勿直接重头运行，可给method赋值为当前报错的算法来debug
  cat(match(method, methods), ":", method, "\n") # 输出当前方法
  method_name = method # 本轮算法名称
  method <- strsplit(method, "\\+")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method)
  
  selected.var = preTrain.var[[method[1]]]
  # 如果筛选出的变量小于阈值，则该算法组合无意义，置空（尤其针对以RSF筛选变量的情况，需在ML脚本中尝试调参）
  if (length(selected.var) <= min.selected.var) {
    model[[method_name]] <- NULL
  } else {
    model[[method_name]] <- RunML(method = method[2], # 用于构建最终模型的机器学习方法
                                  Train_expr = Train_expr[, selected.var], # 训练集有潜在预测价值的变量
                                  Train_surv = Train_surv, # 训练集生存数据
                                  mode = "Model",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                  classVar = classVar)  # 用于训练的生存变量，必须出现在Train_surv中
  }
  
  # 如果最终筛选出的变量小于阈值，则该算法组合也无意义，置空
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
saveRDS(model, file.path(res.path, "model.rds")) # 保存所有模型输出

# 当要求最终模型为多变量cox时，对模型进行更新
if (FinalModel == "multiCox"){
  coxmodel <- lapply(model, function(fit){ # 根据各算法最终获得的变量，构建多变量cox模型，从而以cox回归系数和特征表达计算单样本风险得分
    tmp <- coxph(formula = Surv(Train_surv[[timeVar]], Train_surv[[statusVar]]) ~ .,
                 data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
    return(tmp)
  })
}
saveRDS(coxmodel, file.path(res.path, "coxmodel.rds")) # 保存最终以多变量cox拟合所筛选变量的模型

## Evaluate the model -----------------------------------------------------

# 读取已保存的模型列表（请根据需要调整）
model <- readRDS(file.path(res.path, "model.rds")) # 若希望使用各自模型的线性组合函数计算得分，请运行此行
# model <- readRDS(file.path(res.path, "coxmodel.rds")) # 若希望使用多变量cox模型计算得分，请运行此行

methodsValid <- names(model) # 取出有效的模型（变量数目小于阈值的模型视为无效）

# 根据给定表达量计算样本风险评分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalRiskScore(fit = model[[method]], 
                                    new_data = rbind.data.frame(Train_set,Test_set), # 4.0更新
                                    type = "lp") # 同原文，使用linear Predictor计算得分

}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]]) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){ data.frame(ExtractVar(fit)) }) # 2.1版本更新，提取当B模型依旧降维情况下的最终变量
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"  # 数据框有两列，包含算法以及算法所筛选出的变量
write.table(fea_df, file.path(res.path, "fea_df.txt"),sep = "\t", row.names = F, col.names = T, quote = F)

# 对各模型计算C-index
Cindexlist <- list()
for (method in methodsValid){
  Cindexlist[[method]] <- RunEval(fit = model[[method]], # 预后模型
                                  Test_expr = Test_set, # 测试集预后变量，应当包含训练集中所有的变量，否则会报错
                                  Test_surv = Test_surv, # 训练集生存数据，应当包含训练集中所有的变量，否则会报错
                                  Train_expr = Train_set, # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                  Train_surv = Train_surv, # 若需要同时评估训练集，则给出训练集生存数据，否则置NULL
                                  Train_name = "TCGA", # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                  #Train_expr = NULL,
                                  #Train_surv = NULL, 
                                  cohortVar = "Cohort", # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                  timeVar = timeVar, # 用于评估的生存时间，必须出现在Test_surv中；这里是OS.time
                                  statusVar = statusVar) # 用于评估的生存状态，必须出现在Test_surv中；这里是OS
}
Cindex_mat <- do.call(rbind, Cindexlist)
write.table(Cindex_mat, file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

Cindex_mat <- read.table(file.path(res.path, "cindex_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_Cindex <- sort(apply(Cindex_mat, 1, mean), decreasing = T) # 计算每种算法在所有队列中平均C-index，并降序排列
Cindex_mat <- Cindex_mat[names(avg_Cindex), ] # 对C-index矩阵排序
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3)) # 保留三位小数
fea_sel <- fea_list[[rownames(Cindex_mat)[1]]] # 最优模型（即测试集[或者训练集+测试集]C指数均值最大）所筛选的特征

CohortCol <- brewer.pal(n = ncol(Cindex_mat), name = "Paired") # 设置绘图时的队列颜色
names(CohortCol) <- colnames(Cindex_mat)

# 调用简易绘图函数
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = Cindex_mat, # 主矩阵
                    avg_Cindex = avg_Cindex, # 侧边柱状图
                    CohortCol = CohortCol, # 列标签颜色
                    barCol = "steelblue", # 右侧柱状图颜色
                    col = c("#1CB8B2", "#FFFFFF", "#EEB849"), # 热图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "heatmap of cindex.pdf"), width = cellwidth * ncol(Cindex_mat) + 3, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right") # 热图注释均放在右侧
invisible(dev.off())


#####模型比较
gene<-read.table("Intersect_gene.txt",head=F,sep='\t',check.names = F)
dat<-read.table("GSE17357.txt",head=T,sep='\t',check.names = F,row.names = 1)
expr<-dat[intersect(gene$V1,row.names(dat)),]
dim(expr)
length(gene$V1)
write.table(expr,"GSE17357_expr.txt",quote=F,sep='\t')

work.path <- "I:/Project/结直肠癌-机器学习/17.预后模型比较"; setwd(work.path) 
code.path <- file.path(work.path, "Codes") # 存放脚本
data.path <- file.path(work.path, "InputData") # 存在输入数据（需用户修改）
res.path <- file.path(work.path, "Results") # 存放输出结果
fig.path <- file.path(work.path, "Figures") # 存放输出图片

# 如不存在这些路径则创建路径
if (!dir.exists(data.path)) dir.create(data.path)
if (!dir.exists(res.path)) dir.create(res.path)
if (!dir.exists(fig.path)) dir.create(fig.path)
if (!dir.exists(code.path)) dir.create(code.path)

# 加载R包
library(org.Hs.eg.db)
library(survival)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(openxlsx)

# 加载用于模型比较的脚本
source(file.path(code.path, "compare.R"))
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 训练集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Train_surv <- read.table(file.path(data.path, "Training_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_surv), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_surv <- Train_surv[comsam,,drop = F]

## Validation Cohort -----------------------------------------------------------
# 测试集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与训练集保持相同类型）
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

# 测试集生存数据是行为样本，列为结局信息的数据框（请确保生存时间均大于0）
Test_surv <- read.table(file.path(data.path, "Testing_surv.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_surv), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_surv <- Test_surv[comsam,,drop = F]

# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_surv$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_surv$Cohort, centerFlags = T, scaleFlags = T)
# 读取签名
pubSIG <- read.table(file.path(data.path, "public signatures.txt"), header = T)
if (!"Coef" %in% colnames(pubSIG)) pubSIG$Coef <- NA # 若未匹配到“Coef“列，则新建“Coef“列并默认为NA
pubSIG <- split(pubSIG[, c("SYMBOL", "Coef")], pubSIG$Model)
## My Signature ----------------------------------------------------------------
mySIGname = "MRS" # 本研究所定义的签名的名字，用于在图形中显示
myAlgorithm = "Enet[alpha=0.1]" # 本研究所定义的最优算法，即热图最顶部的算法名称

#--------------------------------------------------------#
# 重要内容请注意，在使用自己定义的签名时有如下两个选项   #
# mySIG1：直接使用根据机器学习算法预测出的风险评分文件   #
# mySIG2：使用机器学习算法定义的基因集拟合多变量Cox计算  #
# 请将对应文件放入Results文件夹，并仅运行对应代码域！！！#
#--------------------------------------------------------#
## 计算C指数 -------------------------------------------------------------------
## mySIG1：RS_mat，使用PrognosticML脚本生成的风险评分文件，此评分是机器学习算法通过predict函数直接获取
mySIG <- read.table(file.path(res.path, "RS_mat.txt"), header = T, check.names = F)
mySIG <- setNames(object = mySIG[[myAlgorithm]], nm = rownames(mySIG))

## mySIG2：fea_df，使用PrognosticML脚本生成的特征文件，并提取基因通过多变量Cox再次计算样本风险评分
mySIG <- read.table(file.path(res.path, "fea_df.txt"), header = T) 
mySIG <- mySIG$features[mySIG$algorithm == myAlgorithm] # 提取PrognosticML脚本在最优算法下获取的特征
mySIG <- data.frame("SYMBOL" = mySIG)

## 整合签名 --------------------------------------------------------------------

# 汇总签名集，以列表的形式构成
# 每一个签名以数据框的形式呈现，至少应当有SYMBOL列，Coef列没有可以填NA
signatures <- pubSIG
signatures[[mySIGname]] <- mySIG

model <- list(); cinfo <- list() # 初始化变量
log.file <- file.path(res.path, "makeCox.log") # 在Results文件夹下新建log文件
if (file.exists(log.file)) file.remove(log.file) # 此log文件用于存放在进行多变量cox分析时的警告
log.file <- file(log.file, open = "a")
sink(log.file, append = TRUE, type = "message")
for (i in names(signatures)){
    if (class(signatures[[i]]) == "data.frame"){
        model[[i]] <- makeCox(Features = signatures[[i]]$SYMBOL, # 签名的基因名
                              coefs = signatures[[i]]$Coef,      # 公共签名所提供的基因系数（如未提供也不必修改此行代码）
                              SIGname = i,                       # 当前循环的签名
                              unmatchR = 0.2,                    # 基因名不匹配率，高于该比率将被剔除；低于匹配率但大于0时会报警告，并存入log文件
                              Train_expr = Train_set,            # 用于计算cox系数的训练集表达谱
                              Train_surv = Train_surv,           # 用于计算cox系数的训练集生存信息
                              statusVar = "OS",                  # 用于构建cox模型的生存结局
                              timeVar = "OS.time")               # 用于构建cox模型的生存时间
    }else{
        model[[i]] = signatures[[i]]
    }
    
    cinfo[[i]] <- calCindex(model = model[[i]],                # 训练的cox模型，为有名字的向量
                            name = i,                          # 当前循环的签名
                            Test_expr = Test_set,              # 用于计算c指数的测试集表达谱
                            Test_surv = Test_surv,             # 用于计算c指数的测试集生存信息
                            Train_expr = Train_set,            # 用于计算c指数的训练集表达谱
                            Train_surv = Train_surv,           # 用于计算c指数的训练集生存信息
                            Train_name = "TCGA",               # 指定训练集的名称
                            #Train_expr = NULL,                # 若不需要评估训练集，则取消此行注释，并注释掉上方对应行
                            #Train_surv = NULL,                # 若不需要评估训练集，则取消此行注释，并注释掉上方对应行
                            CohortVar = "Cohort",              # 用于指定测试集所来自的队列
                            metaCohort = TRUE,                 # 指示是否将测试集合并生成MetaCohort
                            statusVar = "OS",                  # 用于计算c指数的生存结局
                            timeVar = "OS.time")               # 用于计算c指数的生存时间
    message("")
}
closeAllConnections()

cinfo <- do.call(rbind, cinfo)
write.table(cinfo[,1:5], file = file.path(res.path,"cinfo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F) # 输出不同签名在所有队列中的c指数统计量
cinfo <- split(cinfo, cinfo$Cohort)

# 绘图 -------------------------------------------------------------------------

CohortCol <- pal_nejm()(6)
names(CohortCol) <- names(cinfo)

# 批量绘制各个队列的森林图
plots <- lapply(cinfo, function(plot.data){
    plot.data$method <- 
        factor(plot.data$method,
               levels = plot.data$method[order(plot.data$C, decreasing = F)])
    
    # compares two concordance indices: the statistical test is a two-sided Student t test for dependent samples.
    C.compare <- plot.data$C[plot.data$method == mySIGname]
    se.compare <- plot.data$se[plot.data$method == mySIGname]
    n.compare <- plot.data$n[plot.data$method == mySIGname]
    RS.compare <- plot.data$RS[plot.data$method == mySIGname][[1]]
    r.combined <- unlist(lapply(plot.data$RS, function(x) cor(x, RS.compare)))
    var.combined <- plot.data$se^2 + se.compare^2 - 2*r.combined*plot.data$se*se.compare
    p <- pt(abs((plot.data$C-C.compare))/(sqrt(var.combined)), n.compare - 1, lower.tail = F) * 2
    plot.data$label <- cut(p, breaks = c(0, 0.05, 0.01, 0.001, 0.0001))
    plot.data$label <- plyr::mapvalues(x = plot.data$label,
                                       from = c("(0,0.0001]", "(0.0001,0.001]", "(0.001,0.01]", "(0.01,0.05]"), 
                                       to = c("****", "***", "**", "*"))
    
    return(ggplot(plot.data, aes(x = method, y = C, fill = Cohort)) +
               geom_errorbar(aes(ymin = C - 1.96 * se, ymax = C + 1.96 * se), width = .1) +
               geom_point(color = CohortCol[unique(plot.data$Cohort)], size = 2.5) +
               geom_text(aes(x = method, y = max(plot.data$C + 1.96 * plot.data$se - 0.05), label = label)) +
               geom_hline(yintercept = 0.6, linetype = "dashed") +
               ggtitle(label = unique(plot.data$Cohort)) +
               coord_flip() + 
               theme_classic() +
               theme(panel.border = element_rect(fill = NA, size = 1),
                     axis.title = element_blank(),axis.text = element_text(color="black",face="bold",size=12),
                     legend.position = "none"))
})

# 森林图合并并保存
plot_grid(plotlist = plots, nrow = 1)
ggsave(file.path(fig.path, "comparison.pdf"), width = 15, height = 10)
#Survival analysis
custom_theme <- function() {
    theme_survminer() %+replace%
        theme(panel.grid = element_blank(),
              axis.text =   element_text(color='black',size=15, face = "bold"),
              axis.title = element_text(color='black',size=20, face = "bold"),
              panel.background = element_rect(fill = "#EDF8E9")
        )
}
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
    original.p <- p
    if(is.ggplot(original.p)) list.plots <- list(original.p)
    else if(is.list(original.p)) list.plots <- original.p
    else stop("Can't handle an object of class ", class (original.p))
    .set_font <- function(font){
        font <- ggpubr:::.parse_font(font)
        ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
    }
    for(i in 1:length(list.plots)){
        p <- list.plots[[i]]
        if(is.ggplot(p)){
            if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
            if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
            if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
            if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
            if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
            if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
            if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
            list.plots[[i]] <- p
        }
    }
    if(is.ggplot(original.p)) list.plots[[1]]
    else list.plots
}
library(survival)
library("survminer")
rt=read.table("GSE17536.txt",header=T,sep="\t");rt$OS.time=rt$OS.time/365
diff=survdiff(Surv(OS.time,OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(OS.time,OS) ~ risk, data = rt)
p<-ggsurvplot(fit, 
              data=rt,
              conf.int=T,conf.int.style='step',size=1.5,
              pval=paste0 ("P = ",pValue),
              pval.size=12,
              legend.title="risk",
              legend.labs=levels(factor(rt[,"risk"])),
              legend = c(0.9, 0.9),
              font.legend=10,
              xlab="Time(years)",
              break.time.by = 2,
              palette = c("#ee7474","#8bc2d1"),
              surv.median.line = "hv",
              ggtheme = custom_theme(),
              cumevents=F,
              risk.table.height=.30)

p$plot<-customize_labels(
    p$plot,
    font.title    = c(16, "bold", "darkblue"),         
    font.subtitle = c(15, "bold.italic", "black"), 
    font.caption  = c(14, "plain", "black"),        
    font.x        = c(25, "bold", "black"),          
    font.y        = c(25, "bold", "black"),      
    font.xtickslab = c(20, "bold", "black"),
    font.ytickslab = c(20, "bold", "black")
)
#单因素cox分析
library(survival)
rt=read.table("Independence.txt",header=T,sep="\t",check.names=F,row.names=1)            #读取输入文件

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	 cox <- coxph(Surv(OS.time, OS) ~ rt[,i], data = rt)
	 coxSummary = summary(cox)
	 coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	 outTab=rbind(outTab,
	              cbind(id=i,
	              HR=coxSummary$conf.int[,"exp(coef)"],
	              HR.95L=coxSummary$conf.int[,"lower .95"],
	              HR.95H=coxSummary$conf.int[,"upper .95"],
	              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
	              )
}
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

######绘制森林图######
#读取输入文件
rt <- read.table("uniCox.txt",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="Unicox_forest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00A087FF",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "#E64B35FF", "#4DBBD5FF")
points(as.numeric(hr), n:1, pch = 16, col = boxcolor, cex=1.3)
axis(1)
dev.off()
#多因素cox分析
library(survival)
rt=read.table("Multi_independence.txt",header=T,sep="\t",check.names=F,row.names=1)

multiCox=coxph(Surv(OS.time, OS) ~ ., data = rt)
multiCoxSum=summary(multiCox)

outTab=data.frame()
outTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

######绘制森林图######
#读取输入文件
rt <- read.table("multiCox.xls",header=T,sep="\t",row.names=1,check.names=F)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#输出图形
pdf(file="Multi_forest.pdf", width = 7,height = 4)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

#绘制森林图左边的临床信息
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

#绘制森林图
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="#00A087FF",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "#E64B35FF", "#4DBBD5FF")
points(as.numeric(hr), n:1, pch = 16, col = boxcolor, cex=1.3)
axis(1)
dev.off()

#C-index 比较
#计算C-index值
library(survival)
library(survcomp)
dat<-read.table("GSE15459.txt",head=T,sep='\t',check.names = F,row.names = 1)
cox1 <- coxph(Surv(OS.time,OS)~Age,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~Gender,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~Stage,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~riskscore,data = dat,x=TRUE,y=TRUE) 
cindex1 <- concordance.index(predict(cox1),surv.time = dat$OS.time,surv.event = dat$OS,method = "noether")
cindex1$c.index
cindex1$se

#GSEA analysis
library(limma)
rt<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
risk<-read.table('TCGA_risk.txt',head=T,sep='\t',check.names = F,row.names = 1)
rt<-rt[,intersect(row.names(risk),colnames(rt))]
modType=c(rep("High",300),rep("Low",300))
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("High","Low")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="mirnaAll.xls",sep="\t",quote=F)

#GSEA analysis
gsym.fc <- read.table("gene_rank.txt",check.names = F, header = T)
library(clusterProfiler)
library(org.Hs.eg.db)
gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)
#head(gsym.fc.id)

#按照foldchange排序
gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]
#head(gsym.fc.id.sorted)

#获得ENTREZID、foldchange列表，做为GSEA的输入
id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID

kk <- gseKEGG(id.fc, organism = "hsa")
write.table(kk,"KEGG_GSEA.txt",quote=F,sep='\t')
library(GseaVis)
#Low risk
terms <-c("hsa04062","hsa03430","hsa03440","hsa04060","hsa04110","hsa03030","hsa04657")
gseaNb(object = kk,
       geneSetID = terms,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       pFill = "white")
#High risk	   
terms <-c("hsa04010","hsa04015","hsa04514","hsa04151","hsa04330","hsa04510","hsa04512")	   
gseaNb(object = kk,
       geneSetID = terms,
       newGsea = T,
       addPval = T,
       rmHt = T,
       pvalX = 0.9,
       pvalY = 0.5,
       pFill = "white")
	   
#免疫微环境
library(limma)
library(scales)
library(ggplot2)
library(ggtext)
library(reshape2)
library(tidyverse)
library(ggpubr)
riskFile="TCGA_risk.txt"      #风险文件
immFile="infiltration_estimation_for_tcga.csv"     #免疫细胞浸润文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
#读取免疫细胞浸润的文件
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, "riskscore"]
immune=immune[sameSample,]
x=as.numeric(risk)
x[x>quantile(x,0.99)]=quantile(x,0.99)
outTab=data.frame()
for(i in colnames(immune)){
    y=as.numeric(immune[,i])
    if(sd(y)<0.001){next}
    corT=cor.test(x, y, method="spearman")
    cor=corT$estimate
    pvalue=corT$p.value
    if(pvalue<0.05&abs(cor)>0.15){
        outTab=rbind(outTab,cbind(immune=i, cor, pvalue))
        #绘制相关性散点图
        outFile=paste0("cor.", i, ".pdf")
        outFile=gsub("/", "_", outFile)
        df1=as.data.frame(cbind(x,y))
        p1=ggplot(df1, aes(x, y)) + 
            xlab("Risk score") + ylab(i)+
            geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
            stat_cor(method = 'spearman', aes(x =x, y =y))
        #相关性图形
        pdf(file=outFile, width=5, height=4.7)
        print(p1)
        dev.off()
    }
}
#输出相关性结果
write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)
library(ggsci)
corResult=read.table("corResult.txt", head=T, sep="\t")
corResult$Software=sapply(strsplit(corResult[,1],"_"), '[', 2)
corResult$Software=factor(corResult$Software,level=as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))
b=corResult[order(corResult$Software),]
b$immune=factor(b$immune,levels=rev(as.character(b$immune)))
colslabels=rep(pal_npg()(length(levels(b$Software))),table(b$Software))  
ggplot(data=b, aes(x=cor, y=immune, color=Software))+
    labs(x="Correlation coefficient",y="Immune cell")+
    geom_point(size=4.1)+scale_color_npg()+
    theme(panel.background=element_rect(fill="white",size=1,color="black"),
          panel.grid=element_line(color="grey75",size=0.5),
          axis.ticks = element_line(size=0.5),axis.text = element_text(face="bold",size=14),axis.title = element_text(face="bold",size=16),legend.text = element_text(face="bold",size=16),
          axis.text.y = ggtext::element_markdown(colour=rev(colslabels)))

#免疫热图
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskscore[risk$riskscore>quantile(risk$riskscore,0.99)]=quantile(risk$riskscore,0.99)

#读取免疫细胞浸润文件
immune=read.csv(immFile, header=T, sep=",", check.names=F, row.names=1)
immune=as.matrix(immune)

#病人风险值和免疫细胞合并
sameSample=intersect(row.names(risk), row.names(immune))
risk=risk[sameSample, c("risk", "riskscore")]
immune=immune[sameSample,]
data=cbind(risk, immune)

#高低风险组免疫差异分析
outTab=data.frame()
sigCell=c("risk","riskscore")
for(i in colnames(data)[3:ncol(data)]){
    if(sd(data[,i])<0.001){next}
    wilcoxTest=wilcox.test(data[,i] ~ data[,"risk"])
    pvalue=wilcoxTest$p.value
    if(wilcoxTest$p.value<0.01){
        outTab=rbind(outTab,cbind(immune=i, pvalue))
        sigCell=c(sigCell, i)
    }
}
write.table(file="immuneCor.txt", outTab, sep="\t", quote=F, row.names=F)

data=data[,sigCell]
data=data[order(data[,"riskscore"]),]
annCol=data[,1:2]
annCol[,"risk"]=factor(annCol[,"risk"], unique(annCol[,"risk"]))
data=t(data[,(3:ncol(data))])
data <- t(scale(t(data))) 
data[data > 2]  <- 2
data[data < -2] <- -2



annRow=sapply(strsplit(rownames(data),"_"), '[', 2)
annRow=as.data.frame(annRow)
row.names(annRow)=row.names(data)
colnames(annRow)=c("Methods")
annRow[,"Methods"]=factor(annRow[,"Methods"], unique(annRow[,"Methods"]))
gapCol=as.vector(cumsum(table(annCol[,"risk"])))
gapRow=as.vector(cumsum(table(annRow[,"Methods"])))

annColors=list()
annColors[["riskscore"]] <- colorRampPalette(heatmap.BlWtRd)(128)
annColors[["risk"]] <-c("High" = "#4DBBD5FF","Low" = "#E64B35FF")
annColors[["Methods"]] <-c("TIMER" = "#E64B35FF","CIBERSORT" = "#4DBBD5FF","CIBERSORT-ABS"="#00A087FF",
"QUANTISEQ"="#3C5488FF","MCPCOUNTER"="#F39B7FFF","XCELL"="#8491B4FF","EPIC"="#91D1C2FF")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
pheatmap(data,
         annotation=annCol,
         annotation_row=annRow,
         annotation_colors = annColors,
         color = colorRampPalette(heatmap.BlWtRd)(128),
         cluster_cols =F,
         cluster_rows =F,
         gaps_row=gapRow,
         gaps_col=gapCol,
         show_colnames=F,
         show_rownames=T,
         fontsize=15)
	   
#Immune_modulator
tcga_expr <- read.table("Immunegene_expr.txt", row.names = 1, header = T, check.names = F, sep='\t',)
tcga_expr<-t(tcga_expr)
library(biomaRt)
library(pheatmap)
library(RColorBrewer)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
expr2<-tcga_expr
sam.order <- sort(expr2["riskscore",],decreasing = F)
imm.expr <-expr2[row.names(tcga_expr[-1,]), names(sam.order)]
cor.res <- NULL
for (i in 1:nrow(imm.expr)) {
    tmp <- cor.test(as.numeric(imm.expr[i,]),as.numeric(sam.order),method = "spearman") # 使用非参相关性分析
    cor.res <- rbind.data.frame(cor.res,data.frame(gene = rownames(imm.expr)[i], rho = tmp$estimate, p = tmp$p.value, stringsAsFactors = F))
}

pos.cor.gene <- cor.res[which(cor.res$rho > 0.1 & cor.res$p < 0.05), c("gene","rho")] # 显著正相关的基因
neg.cor.gene <- cor.res[which(cor.res$rho < -0.1 & cor.res$p < 0.05), c("gene","rho")] # 显著负相关性的基因（若用例文的-0.4cutoff，这里只能找到2个基因，最终图的效果可能不明显）

# 升序排列，这会让热图呈现出例文那样的楔形渐变的效果
pos.cor.gene <- pos.cor.gene[order(pos.cor.gene$rho,decreasing = F),] 
neg.cor.gene <- neg.cor.gene[order(neg.cor.gene$rho,decreasing = F),]
Sinfo <- read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
Sinfo <- Sinfo[intersect(colnames(expr2),rownames(Sinfo)),]
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
red    <- "#EA6767"
green  <- "#70C17A"
blue   <- "#445EAD"
grey   <- "#7A7A7A"
orange <- "#F7C07C"
yellow <- "#CEDC7C"
Sinfo <- Sinfo[intersect(colnames(expr2),rownames(Sinfo)),]
annCol <- data.frame("riskscore" = as.numeric(sam.order[rownames(Sinfo)]),
                     'Status' = Sinfo$Status,
					 "Age" = Sinfo$Age,
					 'Gender' = Sinfo$Gender,
					 'Stage' = Sinfo$Stage,
					 'T' = Sinfo$T,
					 'N' =Sinfo$N,
					 'M' = Sinfo$M,
                     stringsAsFactors = F,row.names = rownames(Sinfo),
                     check.names = F) 
annCol[is.na(annCol)] <- "NA"

annColors <- list()
annColors[["riskscore"]] <- colorRampPalette(heatmap.BlWtRd)(128)
annColors[["Status"]] <- c("Alive" = "lavenderblush","Dead" = "slategray")
annColors[["Age"]] <- c("<=65" = "#E18727E5",">65" = "#CF4E27")
annColors[["Gender"]] <- c("Female" = '#E0864A',"Male" = 'rosybrown')
annColors[["Stage"]] <- c("Stage I" = "cornsilk","Stage II" = "paleturquoise","Stage III" = "goldenrod","Stage IV" = "firebrick")
annColors[["T"]] <- c("T1" = "cornsilk","T2" = "paleturquoise","T3"='goldenrod',"T4"= 'firebrick','NA'='white')
annColors[["N"]] <- c("N0" = "paleturquoise","N1" = "goldenrod",'N2'='firebrick','NA'='white')
annColors[["M"]] <- c("M0" = "paleturquoise","M1" = "firebrick",'NA'='white')


indata <- t(scale(t(imm.expr[c(pos.cor.gene$gene,neg.cor.gene$gene),names(sam.order)]))) # 注意数据要根据TIM-3的表达顺序排列
indata[indata > 2]  <- 2
indata[indata < -2] <- -2

pheatmap(mat = indata,
         scale = "none", fontsize = 16,
         border_color = NA,
         color = colorRampPalette(heatmap.BlWtRd)(128), # 例文配色
         cluster_cols = F, # 列不聚类
         cluster_rows = F, # 行不聚类
         show_rownames = T, # 不显示行名
         show_colnames = F, # 不显示列名
         annotation_col = annCol[names(sam.order),], # 列注释（注意要根据TIM-3的表达顺序排列）
         annotation_colors = annColors, # 注释颜色
         gaps_row = length(pos.cor.gene$gene), # 在正相关基因结束后分割热图
)
#免疫治疗
generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
    in_gct <- data.frame(GeneID=rownames(in_gct),
                         description="na",
                         in_gct, 
                         stringsAsFactors = F,
                         check.names = F)
    cat("#1.2\n", file = gct_file)
    cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
    cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
    for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
    
    cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
    cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
    cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式 (SKCM)
skcm.immunotherapy.logNC <- read.table("skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.info <- read.table("skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# 创建submap需要的数据格式 (TCGA)
tmp <- read.table("merge.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1) # submap不允许出现flat value, 因此最好选取过滤掉低表达的表达谱，这里使用的数据过滤了超过90%样本表达值均<1的基因
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# 提出亚型的样本，顺序排列
samples.High <- rownames(ann[which(ann$risk == "High"),])
samples.Low <- rownames(ann[which(ann$risk == "Low"),])

sam_info <- data.frame("risk"=c(samples.High,samples.Low),row.names = c(samples.High,samples.Low))
sam_info$rank <- rep(c(1,2),times=c(length(samples.High),length(samples.Low))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# 产生输出数据的文件名
gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- tmp[GENELIST,rownames(sam_info)]  
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"
library(pheatmap)
pheatmap(tmp, cellwidth = 80, cellheight = 100,
         cluster_rows = F,cluster_cols = F,fontsize = 20,
         color = heatmap.YlGnPe[5:1],
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)))
#基因组突变  
library(limma)
library(ComplexHeatmap)
library(tidyverse)
library(magrittr)
library(readxl)
library(stringr)
library(MOVICS)
library(RColorBrewer)
library(ggsci)
rt=read.table("COADREAD_mc3_gene_level.txt",sep="\t",header=T,check.names=F)    #读取输入文件
#数据处理，如果一个基因有多行，取均值
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
mut.binary <- rt[rowSums(rt) > 0,]
sinfo<-read.table("Clinical.txt",head=T,sep='\t',check.names = F,row.names = 1)
comsam <- intersect(rownames(sinfo), colnames(mut.binary))
mut.binary<-mut.binary[,comsam]
sinfo<-sinfo[comsam,]	 
kirc.movics <- list("clust.res" = sinfo,
                    "mo.method" = "CRC")
kirc.movics$clust.res$samID <- rownames(kirc.movics$clust.res)
kirc.movics$clust.res$clust <- sapply(kirc.movics$clust.res$risk,
                                      switch,
                                      "High" = 1, 
                                      "Low" = 2
)
mut.movics <- compMut(moic.res     = kirc.movics,
                      mut.matrix   = mut.binary, # binary somatic mutation matrix
                      doWord       = FALSE, # generate table in .docx format
                      doPlot       = FALSE, # draw OncoPrint
                      freq.cutoff  = 0.03, # keep those genes that mutated in at least 5% of samples
                      p.cutoff     = 0.05, # keep those genes with nominal p value < 0.05 to draw OncoPrint
                      p.adj.cutoff = 2, # don't care about FDR
                      innerclust   = T, # perform clustering within each subtype
)

write.table(mut.movics, "SNP mutation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
SNP<-read.table("SNP.txt",head=T,sep='\t',check.names = F,row.names = 1)
SNP<-SNP[,comsam]
identical(row.names(sinfo),colnames(SNP))
annCol <- sinfo
annColors <- list()
annColors[["risk"]] <- c("High" = "#4DBBD5FF","Low" = "#E64B35FF")
annColors[["Status"]] <- c("Alive" = "lavenderblush","Dead" = "slategray")
annColors[["Age"]] <- c("<=65" = "#E18727E5",">65" = "#CF4E27")
annColors[["Gender"]] <- c("Female" = '#E0864A',"Male" = 'rosybrown')
annColors[["Stage"]] <- c("Stage I" = "cornsilk","Stage II" = "paleturquoise","Stage III" = "goldenrod","Stage IV" = "firebrick")
annColors[["T"]] <- c("T1" = "cornsilk","T2" = "paleturquoise","T3"='goldenrod',"T4"= 'firebrick','NA'='white')
annColors[["N"]] <- c("N0" = "paleturquoise","N1" = "goldenrod",'N2'='firebrick','NA'='white')
annColors[["M"]] <- c("M0" = "paleturquoise","M1" = "firebrick",'NA'='white')
col = c("Mutated" = "#386CB0")
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                  gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    Mutated = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                  gp = gpar(fill = col["Mutated"], col = NA))
    })
op1<-oncoPrint(SNP, get_type = function(x) x,
               alter_fun = alter_fun, col = col,remove_empty_columns = TRUE,column_order = rownames(sinfo),
               
               bottom_annotation = HeatmapAnnotation(df = annCol,col = annColors),column_split = sinfo$risk)   
	   
#CNV
cna <- read.table("all_lesions.conf_99.txt",sep = "\t", row.names = 1,header = T,check.names = F)
cna <- cna[1:72, 
           c(1,9:(ncol(cna)-1))] 
colnames(cna)[2:ncol(cna)] <- substr(colnames(cna)[2:ncol(cna)],1,15) #ID匹配
comsam <- intersect(rownames(sinfo), colnames(cna))
sinfo<-sinfo[comsam,]
cna <- cna[,c("Descriptor",comsam)]
cna$label <- paste0(gsub(" ","",cna$Descriptor),"-", substr(rownames(cna),1,3))
cna <- cna[!duplicated(cna$label),]; rownames(cna) <- cna$label; cna <- cna[,setdiff(colnames(cna),"label")]
cna.modified <- cna[1:nrow(cna),2:ncol(cna)]
cna.binary <- cna.modified[,comsam]
cna.binary[cna.binary > 1] <- 1
kirc.movics <- list("clust.res" = sinfo,
                    "mo.method" = "STAD")
kirc.movics$clust.res$samID <- rownames(kirc.movics$clust.res)
kirc.movics$clust.res$clust <- sapply(kirc.movics$clust.res$risk,
                                      switch,
                                      "High" = 1, 
                                      "Low" = 2
)
mut.movics <- compMut(moic.res     = kirc.movics,
                      mut.matrix   = cna.binary, # binary somatic mutation matrix
                      doWord       = FALSE, # generate table in .docx format
                      doPlot       = FALSE, # draw OncoPrint
                      freq.cutoff  = 0.03, # keep those genes that mutated in at least 5% of samples
                      p.cutoff     = 0.05, # keep those genes with nominal p value < 0.05 to draw OncoPrint
                      p.adj.cutoff = 2, # don't care about FDR
                      innerclust   = T, # perform clustering within each subtype
)

write.table(mut.movics, "CNV mutation.txt",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(cna,"CNV.txt",quote=F,sep='\t')
CNV<-read.table("CNV.txt",head=T,sep='\t',check.names = F,row.names = 1)
CNV<-CNV[,comsam]
identical(row.names(sinfo),colnames(CNV))

annCol <- sinfo
annColors <- list()
annColors[["risk"]] <- c("High" = "#4DBBD5FF","Low" = "#E64B35FF")
annColors[["Status"]] <- c("Alive" = "lavenderblush","Dead" = "slategray")
annColors[["Age"]] <- c("<=65" = "#E18727E5",">65" = "#CF4E27")
annColors[["Gender"]] <- c("Female" = '#E0864A',"Male" = 'rosybrown')
annColors[["Stage"]] <- c("Stage I" = "cornsilk","Stage II" = "paleturquoise","Stage III" = "goldenrod","Stage IV" = "firebrick")
annColors[["T"]] <- c("T1" = "cornsilk","T2" = "paleturquoise","T3"='goldenrod',"T4"= 'firebrick','NA'='white')
annColors[["N"]] <- c("N0" = "paleturquoise","N1" = "goldenrod",'N2'='firebrick','NA'='white')
annColors[["M"]] <- c("M0" = "paleturquoise","M1" = "firebrick",'NA'='white')
col = c("Gain" = pal_nejm()(9)[1],"Loss"=pal_nejm()(9)[2])
alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                  gp = gpar(fill = "#CCCCCC", col = NA))
    },
    # big blue
    Gain = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                  gp = gpar(fill = col["Gain"], col = NA))},
	Loss = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
                  gp = gpar(fill = col["Loss"], col = NA))
    })
op1<-oncoPrint(CNV, get_type = function(x) x,
               alter_fun = alter_fun, col = col,remove_empty_columns = TRUE,column_order = rownames(sinfo),
			   bottom_annotation = HeatmapAnnotation(df = annCol,col = annColors),column_split = sinfo$risk)	   	   
#Drug identification
library(tidyverse) # 用于读取MAF文件
library(ISOpureR) # 用于纯化表达谱
library(impute) # 用于KNN填补药敏数据
library(pRRophetic) # 用于药敏预测
library(SimDesign) # 用于禁止药敏预测过程输出的信息
library(ggplot2) # 绘图
library(cowplot) # 合并图像

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  

dat<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
Sinfo<-read.table("TCGA_risk.txt",head=T,sep='\t',check.names = F,row.names = 1)
sample<-intersect(row.names(Sinfo),colnames(dat))
dat<-dat[,sample]
tumoexpr =dat[rowMeans(dat)>0.5,]
runpure <- F # 如果想运行就把这个改为T
if(runpure) {
    set.seed(123)
    # Run ISOpureR Step 1 - Cancer Profile Estimation
    ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
    # For reproducible results, set the random seed
    set.seed(456);
    # Run ISOpureR Step 2 - Patient Profile Estimation
    ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
    pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}

if(!runpure) {
    pure.tumoexpr <- tumoexpr
}
#加载CTRP细胞系
ctrp.ccl.anno <- read.table("CTRP_ccl_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
ctrp.cpd.anno <- read.table("CTRP_cpd_anno.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T,quote="")
# 2.加载药敏AUC矩阵并进行数据处理
ctrp.auc <- read.table("CTRP_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
prism.auc <- read.delim("PRISM_AUC.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 数据来自https://depmap.org/portal/download/ Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4
prism.ccl.anno <- prism.auc[,1:5] # 前5列为细胞系注释信息
prism.auc <- prism.auc[,-c(1:5)]

## a. 移除缺失值大于20%的药物
ctrp.auc <- ctrp.auc[,apply(ctrp.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(ctrp.auc)]
prism.auc <- prism.auc[,apply(prism.auc,2,function(x) sum(is.na(x))) < 0.2*nrow(prism.auc)]

## b. 移除CTRP数据里源自haematopoietic_and_lymphoid_tissue的细胞系
rmccl <- paste0("CCL",na.omit(ctrp.ccl.anno[,"master_ccl_id"]))
rownames(ctrp.auc) <- paste0("CCL",rownames(ctrp.auc))
ctrp.auc <- ctrp.auc[intersect(rownames(ctrp.auc),rmccl),]

## c. KNN填补缺失值
ctrp.auc.knn <- impute.knn(as.matrix(ctrp.auc))$data

prism.auc.knn <- impute.knn(as.matrix(prism.auc))$data

## d. 数据量级修正（与作者沟通得知）
ctrp.auc.knn <- ctrp.auc.knn/ceiling(max(ctrp.auc.knn)) # 参考Expression Levels of Therapeutic Targets as Indicators of Sensitivity to Targeted Therapeutics (2019, Molecular Cancer Therapeutics)
prism.auc.knn <- prism.auc.knn/ceiling(max(prism.auc.knn))

#药敏预测
# 加载CCLE细胞系的表达谱，作为训练集
ccl.expr <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) 

# 加载基因注释文件，用于基因ID转换
Ginfo <- read.table("overlapTable27.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T) # 参考FigureYa34count2FPKMv2制作的基因注释文件

# 把基因的ensembl ID转换为gene symbol
ccl.expr <- ccl.expr[,-1]; rownames(ccl.expr) <- sapply(strsplit(rownames(ccl.expr),".",fixed = T),"[",1)
comgene <- intersect(rownames(ccl.expr),rownames(Ginfo))
ccl.expr <- ccl.expr[comgene,]
ccl.expr$gene <- Ginfo[comgene,"genename"]; ccl.expr <- ccl.expr[!duplicated(ccl.expr$gene),]; rownames(ccl.expr) <- ccl.expr$gene; ccl.expr <- ccl.expr[,-ncol(ccl.expr)]
keepgene <- apply(ccl.expr, 1, mad) > 0.5 # 保留表达值有效的基因
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1) # 重置细胞系名
trainPtype <- as.data.frame(ctrp.auc.knn)
ccl.name <- ccl.miss <- c() # 替换细胞系名
for (i in rownames(trainPtype)) {
    if(!is.element(gsub("CCL","",i),ctrp.ccl.anno$master_ccl_id)) {
        cat(i,"\n")
        ccl.miss <- c(ccl.miss, i) # 没有匹配到的细胞系
        ccl.name <- c(ccl.name, i) # 插入未匹配的细胞系
    } else {
        ccl.name <- c(ccl.name,  ctrp.ccl.anno[which(ctrp.ccl.anno$master_ccl_id == gsub("CCL","",i)),"ccl_name"]) # 插入匹配的细胞系
    }
}

cpd.name <- cpd.miss <- c() # 替换药物名
for (i in colnames(trainPtype)) {
    if(!is.element(i,ctrp.cpd.anno$master_cpd_id)) {
        cat(i,"\n")
        cpd.miss <- c(cpd.miss, i) # 没有匹配到的药物
        cpd.name <- c(cpd.name, i) # 插入未匹配的药物
    } else {
        cpd.name <- c(cpd.name,  ctrp.cpd.anno[which(ctrp.cpd.anno$master_cpd_id == i),"cpd_name"]) # 插入匹配的药物
    }
}

rownames(trainPtype) <- ccl.name
trainPtype <- trainPtype[setdiff(rownames(trainPtype),ccl.miss),] # 去除未匹配的细胞系
colnames(trainPtype) <- cpd.name
trainPtype <- trainPtype[,setdiff(colnames(trainPtype),cpd.miss)] # 去除未匹配的药物
comccl <- intersect(rownames(trainPtype),colnames(trainExpr)) # 提取有表达且有药敏的细胞系
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5 # 纯化的测试集取表达稳定的基因
testExpr <- pure.tumoexpr[keepgene,]
# 取训练集和测试集共有的基因
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])
outTab <- NULL
# 循环很慢，请耐心
for (i in 1:ncol(trainPtype)) { 
    display.progress(index = i,totalN = ncol(trainPtype))
    d <- colnames(trainPtype)[i]
    tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于CTRP的AUC可能有0值，因此加一个较小的数值防止报错
    
    # 岭回归预测药物敏感性
    ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                    trainingPtype = tmp,
                                    testExprData = testExpr,
                                    powerTransformPhenotype = F,
                                    selection = 1))
    ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
    outTab <- rbind.data.frame(outTab,ptypeOut)
}

dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
keepgene <- apply(ccl.expr, 1, mad) > 0.5
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = T),"[",1)
trainPtype <- as.data.frame(prism.auc.knn)
rownames(trainPtype) <- prism.ccl.anno[rownames(trainPtype),"cell_line_display_name"]
#colnames(trainPtype) <- sapply(strsplit(colnames(trainPtype)," (",fixed = T), "[",1)
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]

# 测试集
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
testExpr <- pure.tumoexpr[keepgene,]
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])


outTab <- NULL

for (i in 1:ncol(trainPtype)) { 
    display.progress(index = i,totalN = ncol(trainPtype))
    d <- colnames(trainPtype)[i]
    tmp <- log2(as.vector(trainPtype[,d]) + 0.00001) # 由于PRISM的AUC可能有0值，因此加一个较小的数值防止报错
    ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                    trainingPtype = tmp,
                                    testExprData = testExpr,
                                    powerTransformPhenotype = F,
                                    selection = 1))
    ptypeOut <- 2^ptypeOut - 0.00001 # 反对数
    outTab <- rbind.data.frame(outTab,ptypeOut)
}

dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
prism.pred.auc <- outTab
top.pps <- Sinfo[Sinfo$riskscore >= quantile(Sinfo$riskscore,probs = seq(0,1,0.1))[10],] # 定义上十分位的样本
bot.pps <- Sinfo[Sinfo$riskscore <= quantile(Sinfo$riskscore,probs = seq(0,1,0.1))[2],] # 定义下十分位的样本
darkblue <- "#4DBBD5FF"
lightblue <- "#E64B35FF"
ctrp.log2fc <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
    display.progress(index = i,totalN = nrow(ctrp.pred.auc))
    d <- rownames(ctrp.pred.auc)[i]
    a <- mean(as.numeric(ctrp.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
    b <- mean(as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
    fc <- b/a
    log2fc <- log2(fc); names(log2fc) <- d
    ctrp.log2fc <- c(ctrp.log2fc,log2fc)
}

candidate.ctrp <- ctrp.log2fc[ctrp.log2fc > 0.1]
prism.log2fc <- c()
for (i in 1:nrow(prism.pred.auc)) {
    display.progress(index = i,totalN = nrow(prism.pred.auc))
    d <- rownames(prism.pred.auc)[i]
    a <- mean(as.numeric(prism.pred.auc[d,rownames(top.pps)])) # 上十分位数的AUC均值
    b <- mean(as.numeric(prism.pred.auc[d,rownames(bot.pps)])) # 下十分位数的AUC均值
    fc <- b/a
    log2fc <- log2(fc); names(log2fc) <- d
    prism.log2fc <- c(prism.log2fc,log2fc)
}
candidate.prism <- prism.log2fc[prism.log2fc > 0.1]

ctrp.cor <- ctrp.cor.p <- c()
for (i in 1:nrow(ctrp.pred.auc)) {
    display.progress(index = i,totalN = nrow(ctrp.pred.auc))
    d <- rownames(ctrp.pred.auc)[i]
    a <- as.numeric(ctrp.pred.auc[d,rownames(Sinfo)]) 
    b <- as.numeric(Sinfo$riskscore)
    r <- cor.test(a,b,method = "pearson")$estimate; names(r) <- d
    p <- cor.test(a,b,method = "pearson")$p.value; names(p) <- d
    ctrp.cor <- c(ctrp.cor,r)
    ctrp.cor.p <- c(ctrp.cor.p,p)
}
candidate.ctrp2 <- ctrp.cor[ctrp.cor < -0.3]  # 这里我调整了阈值，控制结果数目
ctrp.candidate <- intersect(names(candidate.ctrp),names(candidate.ctrp2))

prism.cor <- prism.cor.p <- c()
for (i in 1:nrow(prism.pred.auc)) {
    display.progress(index = i,totalN = nrow(prism.pred.auc))
    d <- rownames(prism.pred.auc)[i]
    a <- as.numeric(prism.pred.auc[d,rownames(Sinfo)]) 
    b <- as.numeric(Sinfo$riskscore)
    r <- cor.test(a,b,method = "pearson")$estimate; names(r) <- d
    p <- cor.test(a,b,method = "pearson")$p.value; names(p) <- d
    prism.cor <- c(prism.cor,r)
    prism.cor.p <- c(prism.cor.p,p)
}
candidate.prism2 <- prism.cor[prism.cor < -0.4]  
prism.candidate <- intersect(names(candidate.prism),names(candidate.prism2))
cor.data <- data.frame(drug = ctrp.candidate,
                       r = ctrp.cor[ctrp.candidate],
                       p = -log10(ctrp.cor.p[ctrp.candidate]))
p1 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
    geom_segment(aes(xend=0,yend=drug),linetype = 2) +
    geom_point(aes(size=p),col = darkblue) +
    scale_size_continuous(range =c(2,8)) +
    scale_x_reverse(breaks = c(0, -0.3, -0.5),
                    expand = expansion(mult = c(0.01,.1))) + #左右留空
    theme_classic() +
    labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
    theme(legend.position = "bottom", 
          axis.line.y = element_blank(),axis.title = element_text(color="black",face="bold",size=16),axis.text = element_text(color="black",face="bold",size=15))

cor.data <- data.frame(drug = prism.candidate,
                       r = prism.cor[prism.candidate],
                       p = -log10(prism.cor.p[prism.candidate]))
cor.data$drug <- sapply(strsplit(cor.data$drug," (",fixed = T), "[",1)

p2 <- ggplot(data = cor.data,aes(r,forcats::fct_reorder(drug,r,.desc = T))) +
    geom_segment(aes(xend=0,yend=drug),linetype = 2) +
    geom_point(aes(size=p),col = darkblue) +
    scale_size_continuous(range =c(2,8)) +
    scale_x_reverse(breaks = c(0, -0.3, -0.5),
                    expand = expansion(mult = c(0.01,.1))) + #左右留空
    theme_classic() +
    labs(x = "Correlation coefficient", y = "", size = bquote("-log"[10]~"("~italic(P)~"-value)")) + 
    theme(legend.position = "bottom", 
          axis.line.y = element_blank(),axis.title = element_text(color="black",face="bold",size=16),axis.text = element_text(color="black",face="bold",size=15))

ctrp.boxdata <- NULL
for (d in ctrp.candidate) {
    a <- as.numeric(ctrp.pred.auc[d,rownames(top.pps)]) 
    b <- as.numeric(ctrp.pred.auc[d,rownames(bot.pps)])
    p <- wilcox.test(a,b)$p.value
    s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
    ctrp.boxdata <- rbind.data.frame(ctrp.boxdata,
                                     data.frame(drug = d,
                                                auc = c(a,b),
                                                p = p,
                                                s = s,
                                                group = rep(c("High riskscore","Low riskscore"),c(nrow(top.pps),nrow(bot.pps))),
                                                stringsAsFactors = F),
                                     stringsAsFactors = F)
}
p3 <- ggplot(ctrp.boxdata, aes(drug, auc, fill=group)) + 
    geom_boxplot(aes(col = group),outlier.shape = NA) + 
    # geom_text(aes(drug, y=min(auc) * 1.1, 
    #               label=paste("p=",formatC(p,format = "e",digits = 1))),
    #           data=ctrp.boxdata, 
    #           inherit.aes=F) + 
    geom_text(aes(drug, y=max(auc)), 
              label=ctrp.boxdata$s,
              data=ctrp.boxdata, 
              inherit.aes=F) + 
    scale_fill_manual(values = c(darkblue, lightblue)) + 
    scale_color_manual(values = c(darkblue, lightblue)) + 
    xlab(NULL) + ylab("Estimated AUC value") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 15,face="bold",color="black"),axis.text.y = element_text(size = 15,face="bold",color="black"),axis.title.y = element_text(size = 16,face="bold",color="black"),
          legend.position = "bottom",
          legend.title = element_blank()) 
dat <- ggplot_build(p3)$data[[1]]

p3 <- p3 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)

prism.boxdata <- NULL
for (d in prism.candidate) {
    a <- as.numeric(prism.pred.auc[d,rownames(top.pps)]) 
    b <- as.numeric(prism.pred.auc[d,rownames(bot.pps)])
    p <- wilcox.test(a,b)$p.value
    s <- as.character(cut(p,c(0,0.001,0.01,0.05,1),labels = c("***","**","*","")))
    prism.boxdata <- rbind.data.frame(prism.boxdata,
                                      data.frame(drug = d,
                                                 auc = c(a,b),
                                                 p = p,
                                                 s = s,
                                                 group = rep(c("High riskscore","Low riskscore"),c(nrow(top.pps),nrow(bot.pps))),
                                                 stringsAsFactors = F),
                                      stringsAsFactors = F)
}
prism.boxdata$drug <- sapply(strsplit(prism.boxdata$drug," (",fixed = T), "[",1)

p4 <- ggplot(prism.boxdata, aes(drug, auc, fill=group)) + 
    geom_boxplot(aes(col = group),outlier.shape = NA) + 
    # geom_text(aes(drug, y=min(auc) * 1.1, 
    #               label=paste("p=",formatC(p,format = "e",digits = 1))),
    #           data=prism.boxdata, 
    #           inherit.aes=F) + 
    geom_text(aes(drug, y=max(auc)), 
              label=prism.boxdata$s,
              data=prism.boxdata, 
              inherit.aes=F) + 
    scale_fill_manual(values = c(darkblue, lightblue)) + 
    scale_color_manual(values = c(darkblue, lightblue)) + 
    xlab(NULL) + ylab("Estimated AUC value") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5,vjust = 0.5,size = 15,face="bold",color="black"),axis.text.y = element_text(size = 15,face="bold",color="black"),axis.title.y = element_text(size = 16,face="bold",color="black"),
          legend.position = "bottom",
          legend.title = element_blank())
dat <- ggplot_build(p4)$data[[1]]

p4 <- p4 + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
plot_grid(p1, p3, p2, p4, labels=c("A", "", "B", ""), 
          ncol=2, 
          rel_widths = c(2, 2))

#drug target,药物细胞系分析
library(readr)
library(ggplot2)
library(ggrepel)
library(cowplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
Tinfo <- read.delim("table s10 target.txt",sep = "\t",row.names = NULL,header = T,check.names = F,stringsAsFactors = F)
Cinfo <- read.csv("sample_info.csv", row.names = 1,check.names = F,stringsAsFactors = F,header = T)
expr.ccle <- read.csv("CCLE_expression.csv", row.names = 1, check.names = F, stringsAsFactors = F, header = T)
ceres <- read.csv("CRISPR_gene_effect.csv", row.names = 1,check.names = F, stringsAsFactors = F, header = T)
Cinfo.hcc <- Cinfo[which(Cinfo$primary_disease == "Colon/Colorectal Cancer"),]
expr.ccle<-expr.ccle[intersect(row.names(Cinfo.hcc),row.names(expr.ccle)),]
colnames(expr.ccle) <- sapply(strsplit(colnames(expr.ccle), " (", fixed = T), "[", 1)
expr.ccle <- as.data.frame(t(expr.ccle))
gene<-read.table("gene.txt",head=F,sep='\t',check.names = F)
ccle.set<-expr.ccle[intersect(gene$V1,row.names(expr.ccle)),]
ccle.set<-t(ccle.set)
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)
model <- readRDS("model.rds") 
methodsValid <- names(model)

source("ML.R")
ccle.set = scaleData(data = ccle.set, centerFlags = T, scaleFlags = T)
RS <- CalRiskScore(fit = model[["Ridge"]], 
                                      new_data = ccle.set, 
                                      type = "lp") 
write.table(RS,"CCLE_risk.txt",quote=F,sep='\t')
#提取共有基因，既来自TCGA，又来自细胞系
TCGA<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
comtarget <- intersect(Tinfo$`Target genes`, rownames(TCGA))
comtarget <- intersect(comtarget, rownames(expr.ccle))
NKS<-read.table("TCGA_Risk.txt",head=T,sep='\t',check.names = F,row.names = 1)
corTCGA <- corCERES <- NULL
for (i in comtarget) {
    
    drug <- paste(Tinfo[which(Tinfo$`Target genes` == i), "Agent names"], collapse = " | ") # 确定该靶点所对应的药物
    
    ## TCGA表达和临床样本的PPS相关性分析
    cor <- cor.test(as.numeric(TCGA[i, names(NKS)]),
                    as.numeric(NKS),
                    method = "pearson") # 原文采用spearman相关性
    corTCGA <- rbind.data.frame(corTCGA,
                                data.frame(target = i,
                                           r = cor$estimate,
                                           p = cor$p.value,
                                           drug = drug,
                                           row.names = i,
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
    
    ## CERES矩阵和细胞系的PPS相关性分析
    cor <- cor.test(as.numeric(expr.ccle[i, names(RS)]),
                    as.numeric(RS),
                    method = "pearson")
    corCERES <- rbind.data.frame(corCERES,
                                 data.frame(target = i,
                                            r = cor$estimate,
                                            p = cor$p.value,
                                            drug = drug,
                                            row.names = i,
                                            stringsAsFactors = F),
                                 stringsAsFactors = F)
    
}
write.table(corTCGA, file = "output_correlation_tcga.txt", sep = "\t",row.names = F, col.names = T, quote = F)
write.table(corCERES, file = "output_correlation_ccle.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# 根据原文阈值筛选对应的候选靶点
candidate.TCGA <- corTCGA[which(corTCGA$r > 0.08 & corTCGA$p < 0.05),] # 蛋白谱靶点值与PPS需正相关
candidate.CERES <- corCERES[which(corCERES$r < -0.236 & corCERES$p < 0.05),] # CERES靶点值与PPS需负相关
candidate.target <- intersect(candidate.TCGA$target, candidate.CERES$target) # 匹配上了原文中的3个，也许还存在数据标准化问题导致结果不完全一致
# 把筛选出的药物靶点保存到文件
write.table(candidate.target, "output_candidate.target.txt", quote = F, row.names = F)
# 设置颜色
grey <- "#BFBFBF"
lightred <- "#FDC9B5"
red <- "#E9583A"
lightblue <- "#66C2A4"
blue <- "#006D2C"
corTCGA$color <- "A" # 给颜色做编号，基础为A
corTCGA[which(corTCGA$r > 0.2 & corTCGA$p < 0.05), "color"] <- "B" # 显著的为B
corTCGA[candidate.target,"color"] <- "C" # 感兴趣的为C
corTCGA$size <- "A" # 同理给大小做编号
corTCGA[which(corTCGA$r > 0.2 & corTCGA$p < 0.05), "size"] <- "B"
corTCGA[candidate.target,"size"] <- "C"
corTCGA <- corTCGA[order(corTCGA$color),] # 让想突出的靶点出现在列表的末尾，这样绘图时不会被其他点所遮挡

selecttargets <- corTCGA[candidate.target,]
selecttargets$label <- rownames(selecttargets)
p1 <- ggplot(corTCGA, aes(r, -log10(p))) + 
    geom_point(aes(color = color, size = size)) + 
    scale_color_manual(values = c(grey, lightred, red))+ # 将编号映射到对应颜色
    scale_size_manual(values = c(1,2,3)) +  # 将编号映射到对应大小
    xlab("Spearman's rank correlation coefficient") + 
    ylab("-log10(P-value)") +
    scale_x_continuous(
        breaks = c(-0.6,-0.3,0,0.3,0.6), # x轴的一些修饰
        labels = c(-0.6,-0.3,0,0.3,0.6),
        limits = c(-0.6, 0.6)) + 
    geom_vline(xintercept = 0.2, color="grey70", # 添加垂直相关性阈值线
               linetype = "longdash", lwd = 0.6) + 
    geom_hline(yintercept = -log10(0.05), color = "grey70", # 添加水平P值阈值线
               linetype = "longdash", lwd = 0.6) +
    theme_bw() +
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = 16, color = "black",face="bold"),
          axis.title = element_text(size = 18, color = "black",face="bold"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
p2 <- p1 + 
    geom_text_repel(data = selecttargets,
                    aes(x = r, y = -log10(p), 
                        label = label),
                    colour="black", nudge_x = .15,
                    box.padding = 0.5,
                    nudge_y = 1,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    segment.angle = 20,
                    size = 10, min.segment.length = 0,
                    point.padding = unit(1, "lines"))

p.tcga <- p2

#CCLE
corCERES$color <- "A"
corCERES[which(corCERES$r < -0.2&corCERES$p < 0.05), "color"] <- "B"
corCERES[candidate.target, "color"] <- "C"
corCERES$size <- "A"
corCERES[which(corCERES$r < -0.2&corCERES$p < 0.05), "size"] <- "B"
corCERES[candidate.target, "size"] <- "C"
corCERES <- corCERES[order(corCERES$color),] # 让想突出的靶点出现在列表的末尾，这样绘图时不会被其他点所遮挡

selecttargets <- corCERES[candidate.target,]
selecttargets$label <- rownames(selecttargets)

p1 <- ggplot(corCERES, aes(r, -log10(p))) + 
  geom_point(aes(color = color, size = size)) + 
  scale_color_manual(values = c(grey, lightblue, blue))+
  scale_size_manual(values = c(1,2,3)) + 
  xlab("Spearman's rank correlation coefficient") + 
  ylab("-log10(P-value)") +
  scale_x_continuous(
    breaks = c(-0.8,-0.5,0,0.5,0.8), 
    labels = c(-0.8,-0.5,0,0.5,0.8),
    limits = c(-0.8, 0.8)) + 
  geom_vline(xintercept = -0.23, color = "grey70", 
             linetype = "longdash", lwd = 0.6) + 
  geom_hline(yintercept = -log10(0.05), color = "grey70", 
             linetype = "longdash", lwd = 0.6) +
  theme_bw() +
  theme(axis.ticks = element_line(size = 0.2, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 16, color = "black",face="bold"),
        axis.title = element_text(size = 18, color = "black",face="bold"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
p2 <- p1 + 
    geom_text_repel(data = selecttargets,
                    aes(x = r, y = -log10(p), 
                        label = label),
                    colour="black", nudge_x = .15,
                    box.padding = 0.5,
                    nudge_y = 1,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    segment.angle = 20,
                    size = 10, min.segment.length = 0,
                    point.padding = unit(1, "lines"))
p.ccle <- p2
p.ccle

scatterProteome <- scatterCERES <- list()
for (i in candidate.target) {
  tmp <- data.frame(var = as.numeric(TCGA[i, names(NKS)]), NKS = as.numeric(NKS))
  cor <- cor.test(tmp$var,
                  tmp$NKS,
                  method = "pearson")
  txt <- paste0("r = ", round(cor$estimate,2), "\n", "P < 0.00001 ") # 构建相关性值的文字标签
  scatterProteome[[i]] <- 
    ggplot(tmp, aes(NKS, var)) + 
    geom_ribbon(stat = "smooth", method = "lm", se = TRUE, # 先画置信区间的彩带以免遮挡散点
                fill = alpha(lightred, 0.6)) + 
    geom_smooth(span = 2, method = lm, color = red, fill = NA) + # 绘制回归线
    geom_point(color = red, size = 2) + 
    xlab("RMS scores") + 
    ylab(paste0("RNA Abundance of ", i)) +
    theme_bw() +
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 18, color = "black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    annotate("text", # 添加相关性的文字标签
             x = min(tmp$NKS), 
             y = max(tmp$var), 
             hjust = 0, fontface = 4, 
             label = txt)
  scatterProteome[[i]]
  ggsave(filename = paste0("scatter plot between RMS and TCGA abundance of ", i, ".pdf"), width = 5, height = 5)
  
  tmp <- data.frame(var = as.numeric(expr.ccle[i, names(RS)]), RS = as.numeric(RS),cell=names(RS))
  cor <- cor.test(tmp$var,
                  tmp$RS,
                  method = "pearson")
  txt <- paste0("r = ", round(cor$estimate,2), "\n", "P = ", round(cor$p.value, 4))
  scatterCERES[[i]] <- 
    ggplot(tmp, aes(RS, var)) + 
    geom_ribbon(stat = "smooth", method = "lm", se = TRUE,
                fill = alpha(lightblue, 0.6)) + 
    geom_smooth(span = 2, method = lm, color = blue, fill = NA) +
    geom_point(color = blue, size = 2) + 
    geom_hline(yintercept = 0, color=red, # 添加CERES的0值水平线
               linetype="longdash", lwd = 0.6) +
    xlab("RMS scores") + 
    ylab(paste0("CERES score of ", i)) +
    theme_bw() + geom_text_repel(aes(label = cell), size = 3,color="#006D2C")+
    theme(axis.ticks = element_line(size = 0.2, color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 18, color = "black"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    annotate("text", 
             x = min(tmp$RS), 
             y = max(tmp$var), 
             hjust = 0, fontface = 4, 
             label = txt)
  scatterCERES[[i]]
  ggsave(filename = paste0("scatter plot between RMS and ceres score of ", i, ".pdf"), width = 5, height = 5)
}

#ELN与GDSC药物相关性，药物机制分析
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)
data<-read.table("merge.txt",head=T,sep='\t',check.names = F,row.names = 1)
data<-as.matrix(data)
GDSC2_Expr=readRDS(file='GDSC2_Expr.rds')
GDSC2_Res=readRDS(file = 'GDSC2_Res.rds')
GDSC2_Res=exp(GDSC2_Res) 

#药物敏感性
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = data,
              batchCorrect = 'eb',      #"eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData')
#1.相关性分析
tcga_gsva <-read.table("Drug_sensitivity.txt",head=T,sep='\t',check.names = F,row.names = 1)
tcga_gsva <-read.table("Drug_sensitivity.txt",head=T,sep='\t',check.names = F,row.names = 1)
tcga_gsva <- tcga_gsva[colnames(tcga_expr),]
index <- rownames(tcga_expr)
y <- as.numeric(tcga_expr)
head(y)
colnames <- colnames(tcga_gsva)
data <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(tcga_gsva[,i]),y, method="spearman")
  data[i,2] <- test$estimate                                            
  data[i,3] <- test$p.value
}
names(data) <- c("symbol","correlation","pvalue")
head(data)
write.table(data, "output_cor.txt", sep = "\t", quote = F, row.names = F)
#Correlation 
df<-read.table("Drug_output.txt",head=T,sep='\t',check.names = F,stringsAsFactors = F)
rownames(df)<-df$Drug

df$Drug <- factor(df$Drug, levels = rownames(df))
#棒棒图
library(viridis)
ggplot(df, aes(Correlation, Drug)) + geom_segment(aes(x=0, xend=Correlation,  y=Drug,yend=Drug,))+
    geom_point(shape=21,aes(size=abs(Correlation),fill=Correlation))+theme_classic()+scale_fill_viridis_b()+
    theme_bw()+theme(axis.text = element_text(color="black",face="bold",size=16),axis.title =element_text(color="black",face="bold",size=16) )+ 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+labs(y="",x="Correlation coefficient")
#2.机制分析
ggplot(df, aes(Drug,Pathway))+geom_point(aes(color=Correlation,size=abs(Correlation)))+theme_classic()+labs(x="",y="")+scale_color_viridis_b()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,color="black",face="bold",size=16),axis.text.y = element_text(color="black",face="bold",size=16))



#单细胞分析
rm(list=ls())
options(width=160)
setwd("/share/Projects/hongkun/project/ExtraProject/SC_CRC")

library(Seurat)
library(SingleR)
library(harmony)
library(clustree)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggsci)
library(RColorBrewer)

#-----------------
#1. 创建数据 
#----------------
count <- read.table("GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt.gz",header=T,check.names=F,row.names=1)
cell <- read.table("GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt.gz",header = T,sep='\t',check.names = F)
cell$Index <- gsub("[-]",".",cell$Index)
rownames(cell) <- cell$Index

scRNA1 <- CreateSeuratObject(count, project = "CRC", 
                           min.cells = 3, 
                           min.features = 200, 
                           assay = "RNA",
                           meta.data = cell)
scRNA1$orig.ident <- "GSE132257"
Idents(scRNA1)='GSE132257' 
scRNA1 <- PercentageFeatureSet(scRNA1, pattern = "^MT-", col.name = "percent.mt")
summary(scRNA1$percent.mt)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   2.217   3.775   5.261   6.949  19.993 
scRNA1 <- subset(scRNA1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
#33694  18409

#-----------------
#2. 分析数据 
#----------------

##scRNA1降维聚类
scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst")
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(scRNA1))

plot1 <- DimPlot(scRNA1, reduction = "pca", group.by="Sample")
plot2 <- ElbowPlot(scRNA1, ndims=50, reduction="pca") 
plotc <- plot1+plot2
#ggsave("1.group_pca.pdf",plotc,width=10,height=4)
ggsave("1.group_pca.pdf",plotc,width=10,height=4)
##harmony
alldata.harmony <- RunHarmony(scRNA1, group.by.vars = "Sample", reduction = "pca", 
    dims.use = 1:50, assay.use = "RNA")

#alldata.harmony <- RunUMAP(alldata.harmony, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony")

pc.num=1:30   #选取前30个主成分
alldata.harmony <- FindNeighbors(alldata.harmony, dims = pc.num) 

for (res in c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0)) {
    alldata.harmony <- FindClusters(alldata.harmony, graph.name = "RNA_snn", resolution = res)
}

##clustree 
p <- clustree(alldata.harmony@meta.data, prefix = "RNA_snn_res.")
# ggsave("2.harmony_clustree.pdf", plot = p, width = 20, height = 16)
ggsave("2.harmony_clustree_normal.pdf", plot = p, width = 20, height = 16)

#选择res=0.6
Idents(alldata.harmony) <- alldata.harmony@meta.data$RNA_snn_res.0.6
sel.clust = "RNA_snn_res.0.6"
alldata.harmony[["seurat_clusters"]] <-Idents(alldata.harmony)
##降维
###tsne
alldata.harmony = RunTSNE(alldata.harmony, dims = pc.num,reduction = "harmony", reduction.name = "tsne_harmony")
#group_by_cluster
plot1 = DimPlot(alldata.harmony, reduction = "tsne_harmony", label=T) 
#group_by_sample
plot2 = DimPlot(alldata.harmony, reduction = "tsne_harmony", group.by='Sample')
#combinate
plotc <- plot1+plot2
ggsave("3.harmony_tsne.pdf", plot = plotc, width = 10, height = 5)

plot3 = DimPlot(alldata.harmony, reduction = "tsne_harmony", group.by='Cell_type') +scale_colour_d3("category20") 
ggsave("3.harmony_tsne_label_celltype.pdf", plot = plot3, width = 7, height = 6)

# plotc <- plot1+plot2+plot3
# ggsave("3.harmony_tsne.pdf", plot = plotc, width = 14, height = 4)

##UMAP
alldata.harmony <- RunUMAP(alldata.harmony, dims = pc.num,reduction = "harmony", reduction.name = "umap_harmony")
#group_by_cluster
plot4 = DimPlot(alldata.harmony, reduction = "umap_harmony", label=T) 
#group_by_sample
plot5 = DimPlot(alldata.harmony, reduction = "umap_harmony", group.by='Sample')
#combinate
plotc <- plot4+plot5
ggsave("4.harmony_umap.pdf", plot = plotc, width = 10, height = 5)

plot6 = DimPlot(alldata.harmony, reduction = "umap_harmony", group.by='Cell_type') +scale_colour_d3("category20") 
ggsave("4.harmony_umap_label_celltype.pdf", plot = plot6, width = 7, height = 6)
# plotc <- plot4+plot5+plot6
# ggsave("4.harmony_umap.pdf", plot = plotc, width = 14, height = 4)

#DEG
markers_genes <- FindAllMarkers(alldata.harmony, logfc.threshold = 0.2, test.use = "wilcox", 
    min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
#top20 <- markers_genes %>% group_by(cluster) %>% top_n(-20, p_val_adj)
write.csv(markers_genes, file = "5.res0.6_cluster_markers.csv")

saveRDS(alldata.harmony,file="GSE132257_seurat_pca.rds")

#SingleR鉴定细胞类型
library(Seurat)
library(SingleR)
library(celldex)
library(patchwork)

scRNA <- readRDS("GSE132257_seurat_pca.rds")

ref <- readRDS("ref.BlueprintEncodeData.Rds")

testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$RNA_snn_res.0.6
cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.main, 
                    de.method = "classic", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

###ref.main
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
celltype


table(scRNA@meta.data$RNA_snn_res.0.6)


#ref.main
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$RNA_snn_res.0.6 == celltype$ClusterID[i]),'celltype_singleR'] <- celltype$celltype[i]}  ##添加细胞类型鉴定结果

table(scRNA@meta.data$celltype_singleR)

saveRDS(scRNA, "GSE132257_seurat_celltype.rds")
               
#已知marker 鉴定
features <- c("EPCAM","KRT18","KRT8","CD3E","CD4","CD8A","MS4A1","HLA-DRA","CD14","LYZ","MS4A7","FCGR3A","CD68","AIF1","IGFBP7","COL1A1","DERL3","MZB1","MS4A2","CPA3")
#20 
p <- FeaturePlot(scRNA, reduction = "umap_harmony", dims = 1:2, features = features, ncol = 4, order = T)
ggsave("6.UMAP_featureplot.pdf", plot = p, width = 16, height = 20) 
p <- FeaturePlot(scRNA, reduction = "tsne_harmony", dims = 1:2, features = features, ncol = 4, order = T)
ggsave("7.tSNE_featureplot.pdf", plot = p, width = 16, height = 20) 

p1 <- DotPlot(scRNA, features = features, group.by = "RNA_snn_res.0.6", assay = "RNA")+
        theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
        ylab("Cell type")+
        xlab("Marker gene")+
        coord_flip()+
        scale_color_viridis()
ggsave("8.res0.6_cluster_dotplot_flip.pdf", plot = p1, width = 8, height = 7)

#根据结果注释细胞类型
scRNA$celltype <- recode(scRNA$RNA_snn_res.0.6,
                             "0" = "CD8+ Tcells", 
                             "1" = "CD8+ Tcells", 
                             "2" = "CD8+ Tcells", 
                             "3" = "Plasma cells",
                             "4" = "B cells", 
                             "5" = "CD8+ Tcells",
                             "6" = "Epithelial cells",
                             "7" = "B cells", 
                             "8" = "Plasma cells",
                             "9" = "Epithelial cells",
                             "10" = "Plasma cells",
                             "11" = "Plasma cells",
                             "12" = "CD8+ Tcells",
                             "13" = "CD8+ Tcells",
                             "14" = "Fibroblast cells", 
                             "15" = "Epithelial cells",
                             "16" = "Epithelial cells",
                             "17" = "Mast cells",
                             "18" = "CD8+ Tcells",
                             "19" = "Fibroblast cells",
                             "20" = "Endothelial cells",
                             "21" = "Endothelial cells",
                             "22" = "CD8+ Tcells",
                             "23" = "Endothelial cells",
							 "24" = "unknown")

scRNA$celltype <- factor(scRNA$celltype,levels=c("CD8+ Tcells","Plasma cells","B cells","Epithelial cells","Fibroblast cells","Mast cells","Endothelial cells","unknown"))
p2 <- DotPlot(scRNA, features = features, group.by = "celltype", assay = "RNA")+
        theme(axis.text.x=element_text(angle=45,vjust = 1,hjust=1))+
        ylab("Cell type")+
        xlab("Marker gene")+
        coord_flip()+
        scale_color_viridis()
ggsave("9.celltype_dotplot_flip.pdf", plot = p2, width = 8, height = 7)


plot7 = DimPlot(scRNA, reduction = "tsne_harmony", group.by='celltype') +scale_colour_d3("category20") 
ggsave("10.harmony_tsne_corrected_celltype.pdf", plot = plot7, width = 7, height = 6)

plot8 = DimPlot(scRNA, reduction = "umap_harmony", group.by='celltype') +scale_colour_d3("category20") 
ggsave("11.harmony_umap_corrected_celltype.pdf", plot = plot8, width = 7, height = 6)

plotc <- plot1+plot2+plot7
ggsave("12.harmony_tsne.pdf", plot = plotc, width = 14, height = 4)

plotc <- plot4+plot5+plot8
ggsave("13.harmony_umap.pdf", plot = plotc, width = 14, height = 4)

saveRDS(scRNA,file="GSE132257_seurat_celltype.rds")
#计算riskscore
scdata<-readRDS("GSE132257_seurat_celltype.rds")
gene<-read.table("gene.txt",head=F,sep='\t',check.names = F)
scset<-scdata@assays$RNA@counts[intersect(gene$V1,row.names(scdata)),]
library(openxlsx)
library(seqinr)
library(plyr)
library(survival)
library(randomForestSRC)
library(glmnet)
library(plsRcox)
library(superpc)
library(gbm)
library(mixOmics)
library(survcomp)
library(CoxBoost)
library(survivalsvm)
library(BART)
library(snowfall)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggrepel)
model <- readRDS("model.rds") 
methodsValid <- names(model)
source("ML.R")
scset = scaleData(data = scset, centerFlags = T, scaleFlags = T)
scset<-t(scset)
RS <- CalRiskScore(fit = model[["Enet[alpha=0.1]"]], 
                                      new_data = scset, 
                                      type = "lp") 
scdata$riskscore <- RS
summary(RS)
#正常和肿瘤样本中细胞组成
table(scdata$celltype)
prop.table(table(scdata$celltype))
table(scdata$celltype, scdata$Class)
Cellratio <- prop.table(table(scdata$celltype, scdata$Class), margin = 2)#计算各组样本不同细胞群比例
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
ggplot(Cellratio) + 
    geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
    scale_fill_manual(values=c(pal_npg()(9)))+
    labs(x='Sample',y = 'Ratio')+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text = element_text(color="black",face="bold",size=16),axis.title = element_text(color="black",face="bold",size=18),legend.title = element_blank())
#表达量点图
scdata2 <- subset(scdata, celltype != "unknown")
DotPlot(object = scdata2, features = c("riskscore"), group.by = 'celltype')+theme(axis.text = element_text(color="black",face="bold",size=14),axis.title = element_text(color="black",face="bold",size=16))
DotPlot(object = scdata2, features = c("TERT"), group.by = 'celltype')+theme(axis.text = element_text(color="black",face="bold",size=14),axis.title = element_text(color="black",face="bold",size=16))