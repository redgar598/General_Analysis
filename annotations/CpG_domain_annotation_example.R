TFBS<-read.table("~/GEO_DNAm/TFBS", header=F)
colnames(TFBS)<-c("bin","chr","chromStart","chromEnd","name","score", "expCount","expNums","expScores")
TFBS <- split(TFBS, TFBS$chr)

## 450K CpG Info
CpG$Chromosome_37<-as.factor(CpG$Chromosome_37)
levels(CpG$Chromosome_37)<-c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9,23,24)
CpG<-split(CpG, CpG$Chromosome_37)

#Call which CpGs are in TFBS
TFBS_450KCpGs<-lapply(1:length(CpG),function(chr) do.call(rbind,lapply(1:nrow(TFBS[[chr]]), function(y) {
  probes_in_TFBS<-CpG[[chr]][which(CpG[[chr]]$Coordinate_37>=TFBS[[chr]]$chromStart[y]& 
                                         CpG[[chr]]$Coordinate_37<=TFBS[[chr]]$chromEnd[y]),]
  if(nrow(probes_in_TFBS)>0){probes_in_TFBS$TFBS<-TFBS[[chr]]$name[y]}else{}
  if(nrow(probes_in_TFBS)>0){probes_in_TFBS$region<-"TFBS"}else{}
  probes_in_TFBS})))
TFBS_450KCpGs<-do.call(rbind, TFBS_450KCpGs) # 1626 CpGs in monocyte enhancers
