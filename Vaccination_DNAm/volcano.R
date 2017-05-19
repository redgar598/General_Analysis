library(ggplot2)
library(RColorBrewer)
library(scales)


makeVolcano<-function(pvalue, deltabeta, dB_threshold, pval_threshold, legend_title, xlimit){
  #VOLCANO PLOT
  volcano<-data.frame(Pvalue=pvalue, Delta_Beta=deltabeta)
  
  #Thresholds 
  dB<-dB_threshold #delta beta cutoff
  Pv<-pval_threshold #Pvalue cutoff
  
  sta_delbeta<-deltabeta[which(pvalue<=pval_threshold)] 
  sta_delbeta<-sta_delbeta[abs(sta_delbeta)>=dB]
  
  print(paste("Hypermethylated", length(sta_delbeta[which(sta_delbeta>=dB)]), sep=": "))
  print(paste("Hypomethylated", length(sta_delbeta[which(sta_delbeta<=(-dB))]) , sep=": "))
  
  volcano<-volcano[complete.cases(volcano),]
  
  ## positive delta beta is hypomethylated (code for volcano should be right now, should colors change?)
  color3<-sapply(1:nrow(volcano), function(x) if(volcano$Pvalue[x]<=Pv){
    if(abs(volcano$Delta_Beta[x])>dB){
      if(volcano$Delta_Beta[x]>dB){"Increased Methylation\n(with Potential Biological Impact)"}else{"Decreased Methylation\n (with Potential Biological Impact)"}
    }else{if(volcano$Delta_Beta[x]>0){"Increased Methylation"}else{"Decreased Methylation"}}}else{"Not Significantly Different"})
  
  volcano$Interesting_CpG3<-color3
  
  
  # COLORS! define here so they are consistent between plots
  # so even if you don't have CpGs in a color catagory the pattern will be maintained
  myColors <- c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue", "grey")
  
  color_possibilities<-c("Decreased Methylation",
                         "Decreased Methylation\n (with Potential Biological Impact)",
                         "Increased Methylation",
                         "Increased Methylation\n(with Potential Biological Impact)",
                         "Not Significantly Different")
  
  names(myColors) <- color_possibilities
  colscale <- scale_color_manual(name = legend_title,
                                 values = myColors, drop = FALSE)
  
  
  #omg
  volcano_plot<-ggplot(volcano, aes(Delta_Beta, -log10(Pvalue), color=Interesting_CpG3))+
    geom_point(shape=19, size=1)+theme_bw()+
    colscale+
    geom_vline(xintercept=c(-dB,dB), color="grey60")+
    geom_hline(yintercept=-log10(Pv), color="grey60")+
    ylab("Multiple Test Corrected P Value (-log10)")+xlab("Delta Beta")+xlim(-xlimit, xlimit)+
    theme(axis.text = element_text(size =14, color="black"),
          axis.title = element_text(size =20),
          legend.text = element_text(size =14),
          legend.title = element_text(size =20))+ 
    guides(color = guide_legend(override.aes = list(size = 4)))
  volcano_plot
}