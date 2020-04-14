library(ggplot2)

options(stringsAsFactors = F)

hikes<-read.csv("hike_plot/hikes.csv")

           # make neater
            hikes<-hikes[which(!(hikes$Name%in%c("Velodrome Trail","Jug Island Beach","High Note Trail","BCMC Trail","UBC Malcolm Knapp Research Forest","Petgill Lake"))),]

rename<-which(!(levels(hikes$Name)%in%c("Velodrome Trail","Jug Island Beach","High Note Trail",
                                        "BCMC Trail","UBC Malcolm Knapp Research Forest","Petgill Lake")))

            #levels(hikes$Name)<-c(levels(hikes$Name)[1:4],"Blackcomb",levels(hikes$Name)[6:9],"Bro. Crk.",
                                   # levels(hikes$Name)[11:22],"Hollyburn",levels(hikes$Name)[24:31],"Seymour",levels(hikes$Name)[33:46])


hikes<-rbind(hikes,list("fimmvorduhals_total","Iceland",
                        "Difficult",12,26,"July – September",0,1300,
                        0, NA, 0))


hikes$steepness<-hikes$Elevation/hikes$Distance


# by steepness
ggplot(hikes, aes(Distance, Elevation)) + 
geom_point(aes(size=Time,
               colour=factor(Completion), 
               fill = steepness), shape=21) +
  geom_text(aes(label=Name),hjust=-0.1, vjust=0.1, color="grey40", size=2.5)+
  scale_color_manual(values=c("white", "red"), name="Completion") + 
  scale_size(range = c(2, 15), name="Time (hours)")+
  scale_fill_gradient(low="#c6dbef", high="#08519c", name="Steepness (m/km)",trans = "sqrt")+
  xlab("Distance (km)")+ylab("Elevation (m)")+theme_bw()

ggsave(file="~/Documents/hikes.pdf")

# by steepness
ggplot(subset(hikes, , aes(Distance, Elevation)) + 
  geom_point(aes(size=Time,
                 colour=factor(Completion), 
                 fill = steepness), shape=21) +
  geom_text(aes(label=Name),hjust=-0.1, vjust=0.1, color="grey40", size=2.5)+
  scale_color_manual(values=c("white", "red"), name="Completion") + 
  scale_size(range = c(2, 15), name="Time (hours)")+
  scale_fill_gradient(low="#c6dbef", high="#08519c", name="Steepness (m/km)",trans = "sqrt")+
  xlab("Distance (km)")+ylab("Elevation (m)")+theme_bw()

# by rating
ggplot(hikes, aes(Distance, Elevation)) + 
  geom_point(aes(size=Time,
                 colour=factor(Completion), 
                 fill = rating), shape=21) +
  geom_text(aes(label=Name),hjust=-0.1, vjust=0.1, color="grey40", size=2.5)+
  scale_color_manual(values=c("white", "red"), name="Completion") + 
  scale_size(range = c(2, 15), name="Time (hours)")+
  scale_fill_gradient(low="#c6dbef", high="#08519c", name="Vancouver Trails Rating")+
  xlab("Distance (km)")+ylab("Elevation (m)")+theme_bw()


####### WCT
          
library(ggplot2)

hikes<-read.csv("hike_plot/hikes_WCT.csv")

hikes$steepness<-hikes$Elevation/hikes$Distance
# make neater
hikes<-hikes[which(!(hikes$Name%in%c("Velodrome Trail","Jug Island Beach","High Note Trail","BCMC Trail","UBC Malcolm Knapp Research Forest","Petgill Lake"))),]

rename<-which(!(levels(hikes$Name)%in%c("Velodrome Trail","Jug Island Beach","High Note Trail",
                                        "BCMC Trail","UBC Malcolm Knapp Research Forest","Petgill Lake")))

hikes<-rbind(hikes,list("fimmvorduhals_total","Iceland",
                        "Difficult",12,26,"July – September",0,1300,
                        0, NA, 0))


hikes$steepness<-hikes$Elevation/hikes$Distance

ggplot(hikes, aes(Distance, Elevation)) + 
  geom_point(aes(size=Time,
                 colour=factor(Completion), 
                 fill = steepness), shape=21) +
  geom_text(aes(label=Name),hjust=-0.1, vjust=0.1, color="grey40", size=2.5)+
  scale_color_manual(values=c("white", "red"), name="Completion") + 
  scale_size(range = c(2, 15), name="Time (hours)")+
  scale_fill_gradient(low="#c6dbef", high="#08519c", name="Steepness (m/km)",trans = "sqrt")+
  xlab("Distance (km)")+ylab("Elevation (m)")+theme_bw()


ggsave(file="~/Documents/hikes_WCT.pdf")
          
