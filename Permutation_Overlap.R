############
## Gene list overlap significant testing
############
# This code takes CpG lists as input, associates the CpGs with genes, 
# asks how much the gene lists overlap (could be easily edited to operate 
# on the CpG level or really any list of things), then permutated 10,000 
# random lists, asks how much then overlap and uses the perumutated overlap
# numbers to ask if the real overlap significantly different than chance overlap



## function is dependent on Rachel's gene annotation
load("Gene_CpG_Relations_updatejune2015.RData")
#https://github.com/redgar598/General_Analysis/tree/master/annotations

# PAWS_Beta is a data frame of beta values

### This function takes two CpG probe lists and asks how much the genes associated with those probes overlap
Permutate_overlap<-function(probe_list1, probe_list2){
  len1<-nrow(probe_list1)
  len2<-nrow(probe_list2)
  rnd1<-rownames(PAWS_Beta)[sample(1:nrow(PAWS_Beta), len1)]
  rnd2<-rownames(PAWS_Beta)[sample(1:nrow(PAWS_Beta), len2)]
  
  Gene1<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rnd1),]
  Gene1<-Gene1[!duplicated(Gene1),]
  Gene1<-Gene1[!duplicated(Gene1[,c(1,4)]),]
  
  Gene2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rnd2),]
  Gene2<-Gene2[!duplicated(Gene2),]
  Gene2<-Gene2[!duplicated(Gene2[,c(1,4)]),]
  
  # if you didn't want to do a gene based overlap you could skip eveything after the first 4 lines of this function and just rest of the function as:
  # length(intersect(rnd1, rnd2)
  length(intersect(unique(Gene1$gene), unique(Gene2$gene))) 
}


# Then apply this function in permutations on random CpG lists
# you put in the the actual CpG lists of interest because the 
# function takes the length of these lists to build the random CpG list
FADV_6_HGHEDLV2_expected_overlap<-sapply(1:100, function(seed){
  set.seed(seed)
  Permutate_overlap(FADV_6_sta_bio_hits, HGHEDLV2_sta_bio_hits)
})

# overlap of real gene hits lists between two variables
# FADV_6$gene is a list og genes associated with the FADV_6_sta_bio_hits CpGs
# or just length(intersect(FADV_6_sta_bio_hits, HGHEDLV2_sta_bio_hits)) 
#if you don't want to look at the gene level
length(intersect(unique(FADV_6$gene), unique(HGHEDLV2$gene))) 

mean(FADV_6_HGHEDLV2_expected_overlap)# 18
sd(FADV_6_HGHEDLV2_expected_overlap)



# fisher's exact test
FADV_6HGHEDLV2_Percent<-barplot$Overlap[1]
FADV_6HGHEDLV2_rnd_Percent<-barplot$Overlap[4]

data<-matrix(c(FADV_6HGHEDLV2_Percent, 100-FADV_6HGHEDLV2_Percent,
               FADV_6HGHEDLV2_rnd_Percent, 100-FADV_6HGHEDLV2_rnd_Percent),
             ncol=2, byrow=T)
FADV_6HGHEDLV2<-fisher.test(round(data))

# Permutation P value
# count the number of permutated random gene lists which overlpa by more on less than your real data
# the code is broken down here for understanding but duplicated in a function aswell
real_overlap<-length(intersect(unique(FADV_6$gene), unique(HGHEDLV2$gene)))
length(which(FADV_6_HGHEDLV2_expected_overlap>real_overlap))/length(FADV_6_HGHEDLV2_expected_overlap)

real_overlap<-length(intersect(unique(deinc2dep$gene), unique(HGHEDLV2$gene))) 
length(which(deinc2dep_HGHEDLV2_expected_overlap>real_overlap))/length(deinc2dep_HGHEDLV2_expected_overlap)

real_overlap<-length(intersect(unique(deinc2dep$gene), unique(FADV_6$gene))) 
length(which(FADV_6_deinc2dep_expected_overlap>real_overlap))/length(FADV_6_deinc2dep_expected_overlap)



### function giving a P value with some interpertation
Overlap_pvalue_function<-function(genelist1, genelist2, permutated_overlaps, multiple_test_correction_number){
  real_overlap<-length(intersect(unique(genelist1), unique(genelist2))) 
  print(paste("Corrected P value for the question are the lists more overlapping than by chance?",
              p.adjust(length(which(permutated_overlaps>=real_overlap))/length(permutated_overlaps), method="fdr", n=multiple_test_correction_number), sep=" "))
  print(paste("Corrected P value for the question are the lists more distinct than by chance?",
              p.adjust(length(which(permutated_overlaps<=real_overlap))/length(permutated_overlaps), method="fdr", n=multiple_test_correction_number), sep=" "))}

## I had three gene lists to compare so a multiple test correction of 3
Overlap_pvalue_function(unique(deinc2dep$gene),unique(FADV_6$gene), FADV_6_deinc2dep_expected_overlap,3)
Overlap_pvalue_function(unique(FADV_6$gene),unique(HGHEDLV2$gene), FADV_6_HGHEDLV2_expected_overlap,3)
Overlap_pvalue_function(unique(deinc2dep$gene),unique(HGHEDLV2$gene), deinc2dep_HGHEDLV2_expected_overlap,3)



