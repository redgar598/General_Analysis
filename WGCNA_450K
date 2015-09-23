setwd("~/WGCNA/")
library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)


## should use regression to adjust for the covariates in the data before performing WGCNA between groups
## also COMBAT should have already been run


# Load Methylation Data
load("~/WGCNA/cell_adjusted_Mcombat_Bmiq_BLBR_Beta_noreplicates_16_All.RData")
load("~/WGCNA/SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]
combat_BLBR_Beta<-combat_BLBR_Beta[which(!(rownames(combat_BLBR_Beta)%in%rownames(SnpatCpG))),] #423384     64
BLBR_Beta<-combat_BLBR_Beta
BLBR_Beta[1:5,1:5]

#rm ba20221 (it is all NAs, was added for the correlations)
BLBR_Beta<-BLBR_Beta[,which(colnames(BLBR_Beta)!="BA20221")]#423384     63

#Select Just Variable CpGs for WGCNA
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
load("~/WGCNA/Variation_measures.RData")## variable acroos all samples
Variation_measures<-Variation_measures[which(!(rownames(combat_BLBR_Beta)%in%rownames(SnpatCpG))),] #423384      5
quantile<-Variation_measures$quantile
BLBR_Beta<-BLBR_Beta[which(quantile>=0.5),]# 19743     63

# transpose the methylation data for further analysis
datMeth0 = as.data.frame(t(BLBR_Beta))
names(datMeth0) = rownames(BLBR_Beta)
rownames(datMeth0) = names(BLBR_Beta)
datMeth0[1:5,1:5]

#Horvath has a step for checking data for excessive missing values and identification of outlier microarray samples
#This should have already been done in quality control steps of methylation data
datMeth<-datMeth0

datMeth_num<-lapply(1:ncol(datMeth), function(y) as.numeric(datMeth[,y]))
datMeth_num<-as.data.frame(do.call(cbind, datMeth_num))
colnames(datMeth_num)<-colnames(datMeth)
rownames(datMeth_num)<-rownames(datMeth)
datMeth<-datMeth_num


#Loading clinical trait data
traitData<-read.csv("~/WGCNA/SampleInfo.csv")
head(traitData)
traitData$Bld_brn<-sapply(1:nrow(traitData), function(x) if(traitData$TissueType[x]=="PBMC"){"Blood"}else{"Brain"})
# remove columns that hold information we do not need.
allTraits = traitData[,];

# blood brain specific allTraits data fixes
fix_names<-unlist(lapply(allTraits$X, function(x) gsub(" ","", x , fixed=TRUE)))
allTraits$X<-fix_names
colnames(allTraits)[1]<-"Sample_ID"
allTraits$Tissue_num<-as.numeric(as.factor(allTraits$TissueType))
allTraits$Bld_brn<-as.numeric(as.factor(allTraits$Bld_brn))



# Form a data frame analogous to expression data that will hold the clinical traits.
allTraits<-allTraits[which(allTraits$Sample_ID%in%rownames(datMeth)),]
datTraits<-allTraits[match(rownames(datMeth), allTraits$Sample_ID),]
datTraits<-datTraits[,c(3,9,10,11,12,14, 15)]

# Re-cluster samples
sampleTree2 = hclust(dist(datMeth), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
#only numeric columns
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datMeth, datTraits, file = "~/WGCNA/BLBR-dataInput_variable.RData")


#######################################################################################################
#Automatic network construction and module detection

load("~/WGCNA/BLBR-dataInput_variable.RData")

#Choosing the soft-thresholding power: analysis of network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datMeth, powerVector = powers, verbose = 7)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

save(sft, file="BLBRvariable-sft.RData")

#One-step network construction and module detection (too slow not enough memory)
net = blockwiseModules(datMeth, power = 7,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "BloodBrainTOMvariable", 
                       verbose = 3)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "BLBRvariable-networkConstruction-auto.RData")


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#######################################################################################################
      # Block-wise network construction and module detection
      bwnet = blockwiseModules(datMeth, maxBlockSize = 10000,
                               power = 4, TOMType = "unsigned", minModuleSize = 30,
                               reassignThreshold = 0, mergeCutHeight = 0.25,
                               numericLabels = TRUE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "BloodBrainTOM-blockwise",
                               verbose = 3)
      table(bwnet$colors)
      # open a graphics window
      sizeGrWindow(12, 9)
      # Convert labels to colors for plotting
      mergedColors = labels2colors(bwnet$colors)
      # Plot the dendrogram and the module colors underneath
      plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                          "Module colors",
                          dendroLabels = FALSE, hang = 0.03,
                          addGuide = TRUE, guideHang = 0.05)
      moduleLabels = bwnet$colors
      moduleColors = labels2colors(bwnet$colors)
      MEs = bwnet$MEs;
      geneTree = bwnet$dendrograms[[1]];
      save(MEs, moduleLabels, moduleColors, geneTree,bwnet,
           file = "~/WGCNA/BLBR-networkConstruction-bwnet.RData")

#######################################################################################################
#Relating modules to external clinical traits
load("BLBRvariable-networkConstruction-auto.RData")

# Define numbers of genes and samples
nGenes = ncol(datMeth);
nSamples = nrow(datMeth);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datMeth, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## Gene relationship to trait and important modules: Gene Significance and Module Membership

# Define variable tissue containing the tissue column of datTrait
weight = as.data.frame(datTraits$Bld_brn);
names(weight) = "Blood_Brain"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datMeth, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datMeth, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#Intramodular analysis: identifying genes with high GS and MM
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "CpG significance for Tissue Type",
                   main = paste("Module membership vs. CpG significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#######################################################################################################
#Summary output of network analysis results
names(datMeth)[moduleColors=="turquoise"]
load("~/WGCNA/Price_annotation.RData")
annotmeth=annotation
probes = names(datMeth)
probes2annot = match(probes, rownames(annotmeth))
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(CpG = probes,
                       geneSymbol = annotmeth$Closest_TSS_gene_name[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Blood_Brain));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_highvar.csv")

## pull one modules genes (this is goo for ermineJ)
dim(geneInfo[which(geneInfo$moduleColor=="turquoise"),])#6850 CpG associated with turquoise
turquoisegeneInfo <- geneInfo[which(geneInfo$moduleColor=="turquoise"),]
write.csv(turquoisegeneInfo, file = "turquoisegeneInfo_highvar.csv")

dim(geneInfo[which(geneInfo$moduleColor=="blue"),])#598 CpG associated with turquoise
bluegeneInfo <- geneInfo[which(geneInfo$moduleColor=="blue"),]
write.csv(bluegeneInfo, file = "bluegeneInfo_highvar.csv")

dim(geneInfo[which(geneInfo$moduleColor=="brown"),])#262 CpG associated with turquoise
browngeneInfo <- geneInfo[which(geneInfo$moduleColor=="brown"),]
write.csv(browngeneInfo, file = "brownegeneInfo_highvar.csv")

dim(geneInfo[which(geneInfo$moduleColor=="yellow"),])#244 CpG associated with turquoise
yellowgeneInfo <- geneInfo[which(geneInfo$moduleColor=="yellow"),]
write.csv(yellowgeneInfo, file = "yellowgeneInfo_highvar.csv")

#######################################################################################################
#Exporting to Cytoscape
load("~/WGCNA/BLBR-dataInput.RData")
load("~/WGCNA/BLBR-networkConstruction-bwnet.RData")



# Read in the annotation file
load("~/WGCNA/Price_annotation.RData")
annotmeth=annotation
# Select modules
modules = c("turquoise");
# Select module probes
probes = names(datMeth)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annotmeth$Closest_TSS_gene_name[match(modProbes, rownames(annotmeth))];



# Restrict the network to the top genes
nTop = 75
IMConn = softConnectivity(datMeth[, modProbes])
top = (order(-IMConn) <= nTop)
datMod<-datMeth[, modProbes]
datModtop<-datMeth[, which(top==T)]


# Recalculate topological overlap if needed
TOM = vectorTOM(datMeth,datModtop, subtract1 = T, blockSize = ncol(datModtop), power = 4);
save(TOM, file="Vecotrized_TOMvariable.RData")

load("Vecotrized_TOMvariable.RData")

# Select the corresponding Topological Overlap
modTOM = TOM[which(top==T),];
dimnames(modTOM) = list(colnames(datModtop), colnames(datModtop))


# Export the network into edge and node list files Cytoscape can read
setwd("~/WGCNA/")
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edgesvariable-", paste(modules, collapse="-"), ".csv", sep=""),
                               nodeFile = paste("CytoscapeInput-nodesvariable-", paste(modules, collapse="-"), ".csv", sep=""),
                               weighted = TRUE,
                               threshold = 0.5,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

## these files are written by horvath's code above which are then read back in a summarized in a way the works faster with cytoscape
edges<-read.csv("CytoscapeInput-edgesvariable-turquoise.csv", sep="")
nodes<-read.csv("CytoscapeInput-nodesvariable-turquoise.csv", sep="")

##Pull out just the gene names columns (keeping the individual CpG names often means the gene are duplicated and it
#takes awhile for cytoscape to load 3000+ rows)
#these lines minimize to one entry per gene
edges<-edges[,3:6]
nodes<-nodes[,2:3]
edges<-edges[!duplicated(edges), ]
nodes<-nodes[!duplicated(nodes), ]
nodes<-nodes[which(nodes$altName%in%c(edges$fromAltName, edges$toAltName)),]

## number of connections
conn<-tapply(edges$fromAltName, edges$toAltName, length)
conn<-data.frame(gene=names(conn), conn=conn)
conn2<-tapply(edges$toAltName, edges$fromAltName, length)
conn2<-data.frame(gene=names(conn2), conn=conn2)
conn<-rbind(conn,conn2)
conn$conn<-as.numeric(conn$conn)
conn<-tapply(conn$conn,conn$gene, sum, na.rm=T)
conn<-data.frame(gene=names(conn), conn=conn)

nodes<-merge(nodes, conn, by.x="altName", by.y="gene")
write.csv(edges,"~/WGCNA/CytoscapeInput-edges-turquoise-highvar.csv", quote=F)
write.csv(nodes,"~/WGCNA/CytoscapeInput-nodes-turquoise-highvar.csv", quote=F)

