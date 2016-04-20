WGCNA in PAWS
========================================================

Starting from QC filtered, nromalized, combatted, cell type corrected beta value data this code:

*Filters SNP CpGs
*Filters non-variable CpGs
*Regresses out covariates in the data
*Runs WGCNA
*Associates the modules with meta data variables

# Libraries
```{r}
setwd("/big_data/redgar/cake_backup/redgar/WGCNA")
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(gridExtra)
library(RCurl)
library(WGCNA)
enableWGCNAThreads() ## doesnt work through R studio. 
                      #Run the blockwiseModules() step in terminal or basic R with enableWGCNAThreads() if possible.
```

# Load data
```{r}
load("~/PAWS/PAWS_meta_sentrix_genetic_clusters.RData")
load("~/PAWS/combat_PAWS_Beta_norep.RData")
meta<-meta[which(meta$Meth_ID%in%colnames(PAWS_Beta)),]
meta<-meta[match(colnames(PAWS_Beta), meta$Meth_ID),]
```

# Remove SNP CpGs
```{r}
load("~/PAWS/SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
PAWS_Beta<-PAWS_Beta[which(!(rownames(PAWS_Beta)%in%rownames(SnpatCpG))),]
PAWS_Beta<-as.data.frame(PAWS_Beta)
```

# Convert to M values
```{r}
#M value transformation
Mval<-function(beta) log2(beta/(1-beta))
PAWS_Mval = apply(PAWS_Beta, 1, Mval) # need mvalues for combat
PAWS_Mval = as.data.frame(PAWS_Mval)
PAWS_Mval = t(PAWS_Mval)
```

## which CpGs are non-variable in buccal cells
```{r}
x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Buccal_CpGs.csv")
y <- read.csv(text = x)

PAWS_independent_buccal_invariable<-PAWS_Beta[which(rownames(PAWS_Beta)%in%y$CpG),]#114770/120009 of the independnt invariable sites are in PAWS

# Call varibility in PAWS
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
PAWS_ref_range<-sapply(1:nrow(PAWS_independent_buccal_invariable), function(x) Variation(PAWS_independent_buccal_invariable[x,]))
Invariable_in_PAWS<-PAWS_independent_buccal_invariable[which(PAWS_ref_range<0.05),]

# Which CpGs are invariable in PAWS and the independent data
invar_in_PAWS_and_independent<-intersect(y$CpG, rownames(Invariable_in_PAWS)) #114211/114770 (99.5%)
PAWS_Beta_variable<-PAWS_Beta[which(!(rownames(PAWS_Beta)%in%invar_in_PAWS_and_independent)),]#295667 
PAWS_Beta<-PAWS_Beta_variable #409878 vs 295667 sites

PAWS_Mval_variable<-PAWS_Mval[which(!(rownames(PAWS_Mval)%in%invar_in_PAWS_and_independent)),]#295667 
PAWS_Mval<-as.data.frame(PAWS_Mval_variable) #409878 vs 295667 sites
```


## Regress out the co-variates using the residuals
```{r}
## remove samples with missing meta data as can not calculate a residual for them (Using 175 samples for WGCNA)
meta_complete<-meta[complete.cases(meta[,c("Age_Genetic_Collection","minor_child","Genetic_cluster")]),]
  
PAWS_Mval_complete<-PAWS_Mval[,which(colnames(PAWS_Mval)%in%meta_complete$Meth_ID)] # 295667    175
meta_complete<-meta_complete[match(meta_complete$Meth_ID, colnames(PAWS_Mval_complete)),]


# Impute missing mvalues
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
Mval_imputed<-t(imputeMedianv3(PAWS_Mval_complete))

avebeta.lm<-lapply(1:nrow(Mval_imputed), function(x){
  lm(unlist(Mval_imputed[x,])~as.factor(meta_complete$Genetic_cluster)+as.factor(meta_complete$minor_child)+meta_complete$Age_Genetic_Collection)
})

residuals<-lapply(avebeta.lm, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
adj.residuals<-residuals+matrix(apply(Mval_imputed, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

PAWS_covar_adjusted_mval<-adj.residuals
rownames(PAWS_covar_adjusted_mval)<-rownames(Mval_imputed)

save(PAWS_covar_adjusted_mval, file="PAWS_mvalues_adjusted_covariates.RData")
```



# sanity check that there are no associations with the adjusted covariates
```{r}
# ANOVA Genetic Cluster
Genetic_cluster_aov_unadjusted<-sapply(1:nrow(Mval_imputed), function(CpG){
  x<-aov(Mval_imputed[CpG,]~as.factor(meta_complete$Genetic_cluster))
  summary(x)[[1]][["Pr(>F)"]][1]})

Genetic_cluster_aov_adjusted<-sapply(1:nrow(PAWS_covar_adjusted_mval), function(CpG){
  x<-aov(PAWS_covar_adjusted_mval[CpG,]~as.factor(meta_complete$Genetic_cluster))
  summary(x)[[1]][["Pr(>F)"]][1]})


# lm Age
Age_lm_unadjusted<-sapply(1:nrow(Mval_imputed), function(CpG){
  x<-lm(Mval_imputed[CpG,]~meta_complete$Age_Genetic_Collection)
  summary(x)$coefficients[2,"Pr(>|t|)"]})

Age_lm_adjusted<-sapply(1:nrow(PAWS_covar_adjusted_mval), function(CpG){
  x<-lm(PAWS_covar_adjusted_mval[CpG,]~meta_complete$Age_Genetic_Collection)
  summary(x)$coefficients[2,"Pr(>|t|)"]})


# t.test minor child
MC_T_unadjusted<-sapply(1:nrow(Mval_imputed), function(CpG){
  x<-t.test(Mval_imputed[CpG,]~as.factor(meta_complete$minor_child))
  x$p.value})
MC_T_adjusted<-sapply(1:nrow(PAWS_covar_adjusted_mval), function(CpG){
  x<-t.test(PAWS_covar_adjusted_mval[CpG,]~as.factor(meta_complete$minor_child))
  x$p.value})

# t.test gender
gender_T_unadjusted<-sapply(1:nrow(Mval_imputed), function(CpG){
  x<-t.test(Mval_imputed[CpG,]~as.factor(meta_complete$Gender_final))
  x$p.value})
gender_T_adjusted<-sapply(1:nrow(PAWS_covar_adjusted_mval), function(CpG){
  x<-t.test(PAWS_covar_adjusted_mval[CpG,]~as.factor(meta_complete$Gender_final))
  x$p.value})


## Save the p values
save(MC_T_adjusted,MC_T_unadjusted, 
     Age_lm_adjusted, Age_lm_unadjusted, 
     Genetic_cluster_aov_adjusted, Genetic_cluster_aov_unadjusted, 
     gender_T_adjusted, gender_T_unadjusted, file="Pvalues_adjusted_unadjusted.RData")


## plot the p value distributions
pvalue_dist<-data.frame(CpG=rep(c(rownames(PAWS_covar_adjusted_mval),rownames(Mval_imputed)),times=4), 
                        Nominal_P=c(MC_T_adjusted,MC_T_unadjusted, 
                                    Age_lm_adjusted, Age_lm_unadjusted, 
                                    Genetic_cluster_aov_adjusted, Genetic_cluster_aov_unadjusted,
                                    gender_T_adjusted, gender_T_unadjusted),
                        Methylation_values=rep(rep(c("Adjusted","Unadjusted"), each=nrow(Mval_imputed)), times=4),
                        Covariate=rep(c("Minor Child","Age","Genetic Cluster","Gender"), each=nrow(Mval_imputed)*2))

pvalue_dist$Methylation_values<-factor(pvalue_dist$Methylation_values, levels=c("Unadjusted", "Adjusted"))

ggplot(pvalue_dist, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+
  facet_grid(Methylation_values~Covariate, scales="free_y")


## Save the adjusted data and relevant meta data
save(PAWS_covar_adjusted_mval, meta_complete, file="PAWS_covariate_adjusted_Mvalues.RData")
```






# Weighted Gene Co-expression (co-methylation) Network Analysis

From the WGCNA online tutorial "Dealing with large data sets: block-wise network construction and
module detection"

# load and transfrom data
```{r}
load("PAWS_covariate_adjusted_Mvalues.RData")
datMeth<-t(PAWS_covar_adjusted_mval)
```

#Automatic block-wise network construction and module detection


##Choosing the soft-thresholding power: analysis of network topology

Here we define the power to which the correlation matrix values are taken to, to form a adjacency matrix. This power value deflates the importance of weakly correlated CpGs and increases the weight of highly correlated CpGs. The default in WGCNA is a signed network which will also devalue negativly correlated CpGs, and I don't know if this makes sense for methylation. Maybe it does?

But here we are calculating the B in this formula 0.5+0.5 * cor^B . Where cor is the CpG-CpG correlation. The resulting value once the power (B) is chosen is the adjacnecy  value for two CpGs.

This section of code was run in terminal with the multicores lines used. I may do that for all sections if things continue to run slow through Rstudio (can't do the WGCNA multi core line in R studio)

```{r}
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datMeth, powerVector = powers, verbose = 5)
save(sft, file="sft_PAWS.RData")

load("sft_PAWS.RData")

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"), ylim=c(0, 1));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

It seems like a soft threshold of 5 will be appropriate for this analysis

#Block-wise network construction and module detection
```{r}
bwnet = blockwiseModules(datMeth, maxBlockSize = 20000,
                         power = 5, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "PAWS_TOM_blockwise",
                         verbose = 3)
```

#We now save the module assignment and module eigengene information necessary for subsequent analysis.
```{r}
moduleLabels = bwnet$colors
moduleColors = labels2colors(bwnet$colors)
MEs = bwnet$MEs;
geneTree = bwnet$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,bwnet,
file = "PAWS-networkConstruction-blockwise.RData")
```



# Associate modules with traits
```{r}
load("PAWS-networkConstruction-blockwise.RData")
load("~/PAWS/combat_PAWS_Beta_norep.RData")
meta_complete<-meta_complete[which(meta_complete$Meth_ID%in%colnames(PAWS_covar_adjusted_mval)),]
meta_complete<-meta_complete[match(colnames(PAWS_covar_adjusted_mval), meta_complete$Meth_ID),]

datTraits<-meta_complete[c("Genetic_cluster","Gender_final","DECHIETH","Age_Genetic_Collection","COMPSES2","deinc2dep",
                           "HGHEDLV2","BMI_Percentile","Sentrix_ID","FADV_6","FLadder","minor_child")]

# Define numbers of genes and samples
nGenes = ncol(datMeth);
nSamples = nrow(datMeth);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datMeth, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
```


## plot
```{r}
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 11.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = greenWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
```




