### WGCNA Differential heat map
library(ggplot2)


load("~/Documents/Blood_Brain/cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData") #441198     64
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta
load("~/Documents/Blood_Brain/SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[which(!(rownames(combat_BLBR_Beta_adjusted)%in%rownames(SnpatCpG))),]
combat_BLBR_Beta_adjusted<-combat_BLBR_Beta_adjusted[,1:63]


meta<-read.csv("~/Documents/Blood_Brain/SampleInfo.csv")
fix_names<-unlist(lapply(meta$X, function(x) gsub(" ","", x , fixed=TRUE)))
meta$X<-fix_names
meta<-meta[which(meta$X%in%colnames(combat_BLBR_Beta_adjusted)),]
meta<-meta[match(colnames(combat_BLBR_Beta_adjusted), meta$X),]


setwd("~/Documents/Blood_Brain/WGCNA")
turquoise<-read.csv("turquoisegeneInfo_highvar.csv")
load("Eigengenes_highly_variable.RData")

meta_datTraits<-cbind(meta, MEs0)
meta_datTraits<-meta_datTraits[order(meta_datTraits$TissueType),]

meta_datTraits$SubjectNumC<-as.factor(meta_datTraits$SubjectNumC)
meta_datTraits$TissueType<-as.factor(meta_datTraits$TissueType)

ggplot(meta_datTraits, aes(TissueType, MEturquoise, fill=TissueType, group=SubjectNumC))+
  geom_bar(stat="identity", position="dodge", color="black")+theme_bw()+
  scale_fill_manual(values=c("#fb6a4a","#ef3b2c","#cb181d","cornflowerblue"))


############################ heat map with betas of CpGs in module
library("gplots")
library("RColorBrewer")
cols <- colorRampPalette(brewer.pal(10, "BuPu"))(256)

turquoise_top<-turquoise[1:100,]
turquoise_beta<-as.matrix(combat_BLBR_Beta_adjusted[which(rownames(combat_BLBR_Beta_adjusted)%in%turquoise_top$CpG),])
turquoise_beta<-turquoise_beta[,order(datTraits$Tissue_num)]

heatmap.2(turquoise_beta, col=cols, 
          ColSideColors=as.character(datTraits$Tissue_num[order(datTraits$Tissue_num)]),
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexRow=0.5,
          dendrogram="none",labRow=NA,Rowv = TRUE, Colv=FALSE)#
