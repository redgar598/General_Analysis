load("Gene_CpG_Relations_updatejune2015.RData")

## genes associated with hits
FADV_6<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(FADV_6_sta_bio_hits)),]


############# Much more complicated but sooooo worth it gene summary table
## CpG number in a gene adjustment
Gene_CpG_Relations_update$gene<-as.character(Gene_CpG_Relations_update$gene)
Overrep<-as.data.frame(tapply(Gene_CpG_Relations_update$Probe_ID, Gene_CpG_Relations_update$gene, length))
Overrep$Gene<-rownames(Overrep)
colnames(Overrep)<-c("CpG_number", "Gene")
Overrep<-Overrep[which(Overrep$Gene!="None"),]
mean(Overrep$CpG_number, na.rm=T)# 25
Overrep$Enrichment_fromAverage<-Overrep$CpG_number/mean(Overrep$CpG_number, na.rm=T)

## Gene summaries
load("~/BLBR/Price_annotation.RData")
annotation$CpG<-rownames(annotation)

Format_gene_table<-function(Gene_CpG_Relations_update_subset){
  print(paste("CpGs Associated: ", length(unique(Gene_CpG_Relations_update_subset$Probe_ID)), sep=""))
  print(paste("Genes Associated: ", length(unique(Gene_CpG_Relations_update_subset$gene)), sep=""))
  Overrep_subset<-as.data.frame(tapply(Gene_CpG_Relations_update_subset$Probe_ID, Gene_CpG_Relations_update_subset$gene, length))
  Overrep_subset$Gene<-rownames(Overrep_subset)
  colnames(Overrep_subset)<-c("CpG_number", "Gene")
  Overrep_subset<-Overrep_subset[which(Overrep_subset$Gene!="None"),]
  Overrep_subset_merge<-merge(Overrep_subset, Overrep, by="Gene")
  colnames(Overrep_subset_merge)<-c("Gene","CpG_Associated","CpG_in_Gene", "Enrichment_fromAverage")
  Overrep_subset_merge$Suprise<-Overrep_subset_merge$CpG_Associated/Overrep_subset_merge$Enrichment_fromAverage
  Gene_table<-merge(Gene_CpG_Relations_update_subset, Overrep_subset_merge, by.x="gene", by.y="Gene")
  Gene_table<-merge(Gene_table, annotation[,c(49,50,58)], by.x="Probe_ID", by.y="CpG")
  ## these lines add a p value and delta beta to the gene summary table
  ## PAWS_Beta is a a data frame a beta values
  ## if you have an object with multiple test correction p values call it Multi_test_corr_relaxed
  ## if you have an object with delta betas call it delbeta
  ## both Multi_test_corr_relaxed and delbeta need to be the same length and in the same order as PAWS_Beta
  pval<-data.frame(CpG=rownames(PAWS_Beta)[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)],
                   corr_pval=Multi_test_corr_relaxed[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)],
                   db=delbeta[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)])
  Gene_table<-merge(Gene_table, pval, by.x="Probe_ID", by.y="CpG")
  Gene_table<-Gene_table[,c(2,7,8,9,10,1,4,3,6,5,11,12,13,14)]
  Gene_table<-Gene_table[order(-Gene_table$Suprise, Gene_table$gene),]
  Gene_table}




#FADV_6_sta_bio_hits is a dataframe of beta values that are significant with CpGs as rownames (but all you need is a list of CpG names) 
FADV_6<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(FADV_6_sta_bio_hits)),]
FADV_6<-FADV_6[!duplicated(FADV_6),]
FADV_6<-FADV_6[!duplicated(FADV_6[,c(1,4)]),]#remove duplicate CpG to gene associations

FADV_6_genes<-Format_gene_table(FADV_6)