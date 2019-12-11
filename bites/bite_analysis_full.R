library(dplyr)
library(ggplot2)
library(reshape)
library(scales)
library(lubridate)

options(stringsAsFactors = FALSE)
bites <- read.csv("bites/Bite_Tracker - 2017_cleaned.csv")


# sort out some labeling inconsistencies
bites<-as.data.frame(bites)
bites$Date<-gsub("Dec.","Dec",bites$Date)
bites$Date<-as.Date(bites$Date,format='%b %d')
bites$month<-strftime(bites$Date, "%m")
year(bites$Date) <- 2017

bites$Victim<-gsub(" ","",bites$Victim)
bites$Condition<-as.factor(bites$Condition)
levels(bites$Condition)<-c("Drinking","Drinking","Drinking","Drinking","Not Drinking","Not Drinking", "Not Drinking", "Not Drinking")

####################
## PLOTS
####################

#### Seasonality of biting
month_count<-as.data.frame(tapply(bites$Date, bites$month, length))
colnames(month_count)<-"count"
month_count$month<-rownames(month_count)
month_count$month<-as.numeric(month_count$month)
month_count$month_alpha<-c("Jan","Feb","Mar","Apr","May","June","July","Aug","Sep","Oct","Nov","Dec")

ggplot(month_count, aes(reorder(month_alpha,month), count))+
  geom_bar(stat = "identity", col="black",fill="#41ab5d")+theme_bw()+
  xlab("")+ ylab("Total Bite Count")


### Victim
person_count<-as.data.frame(tapply(bites$Date, bites$Victim, length))
colnames(person_count)<-"count"
person_count$victim<-rownames(person_count)

ggplot(person_count, aes(victim, count, fill=victim))+
  geom_bar(stat = "identity", color="black")+theme_bw()+
  scale_fill_manual(values=c("#a6d96a","cornflowerblue"), guide=F)+
  theme(text = element_text(size=20)) +geom_text(aes(label=count), vjust=-0.5, size=6)+
  ylab("Total Bite Count")+xlab("Victim")+ylim(0,45)

ggsave(file="bites/figs/Victim.jpeg", width=5, height=5)

### Drinking overall
drink_count<-as.data.frame(tapply(bites$Date, bites$Condition, length))
colnames(drink_count)<-"count"
drink_count$condition<-rownames(drink_count)

ggplot(drink_count, aes(condition, count, fill=condition))+
  geom_bar(stat = "identity", color="black")+theme_bw()+
  scale_fill_manual(values=c("#fe9929","#d9d9d9"), guide=F)+
  theme(text = element_text(size=20)) +geom_text(aes(label=count), vjust=-1, size=6)+
  ylim(0,40)+ylab("Total Bite Count")+xlab("Condition")

ggsave(file="bites/figs/Drinking.jpeg", width=5, height=5)

### Drinking and Person
drinkvictim_count<-as.data.frame(tapply(bites$Date, list(bites$Condition, bites$Victim), length))
drinkvictim_count$condition<-rownames(drinkvictim_count)

drinkvictim_count<-melt(drinkvictim_count)
drinkvictim_count$col<-paste(drinkvictim_count$condition,drinkvictim_count$variable, sep="-")

ggplot(drinkvictim_count, aes(variable, value, fill=col))+
  geom_bar(stat = "identity", color="black",position = "dodge")+theme_bw()+
  scale_fill_manual(values=c("#238b45","#2171b5","#a1d99b","#9ecae1"), name="")+
  theme(text = element_text(size=20)) +geom_text(aes(label=value), vjust=-1, size=6, position = position_dodge(0.9))+
  ylim(0,28)+ylab("Total Bite Count")+xlab("Victim")

ggsave(file="bites/figs/DrinkingVictim.jpeg", width=10, height=5)


###

### Biting hour
bites$Time2<-hour(hm(bites$Time))
bites$Time2<-as.numeric(gsub(":","",bites$Time))


ggplot(bites, aes(Time2))+
  geom_density(fill="lightgrey", size=1, adjust=1/4)+theme_bw()+
  scale_x_continuous(breaks = seq(0,2400, 300), labels = c("0:00","3:00","6:00","9:00","12:00","15:00","18:00","21:00","24:00"), name="Time")

ggsave(file="bites/figs/BitingHour.jpeg", width=10, height=3)


### git commit style

dates_all<-data.frame(bite=0,Date=seq(min(bites$Date), max(bites$Date), 1))
dates_all$bite<-sapply(1:nrow(dates_all), function(x) if(dates_all$Date[x]%in%bites$Date){length(which(bites$Date==dates_all$Date[x]))}else{0})


dates_all$week<-strftime(dates_all$Date, format = "%V")
dates_all$day<- weekdays(as.Date(dates_all$Date))

dates_all$day<-as.factor(dates_all$day)

#dates_all$day<-factor(dates_all$day,levels(dates_all$day)[c(2,6,7,5,1,3,4)])
dates_all$day<-factor(dates_all$day,levels(dates_all$day)[c(4,3,1,5,7,6,2)])

#month label
dates_all$month<-strftime(dates_all$Date, "%m")
monthbreaks<-as.numeric(tapply(dates_all$week, dates_all$month, min))
dates_all$week<-as.numeric(dates_all$week)
# 
ggplot(dates_all, aes(week,day, fill = as.factor(bite))) +
  geom_tile(color = "black",size=0.5) +
  theme_gray(8)+scale_fill_manual(values=c("#d9f0a3","#41ab5d","#005a32"), name="Bite Count")+
  theme(axis.text = element_text(size =10, color="black"),
        axis.text.x = element_text(),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        axis.ticks = element_blank())+
  xlab("")+ylab(NULL)+
  scale_x_continuous(
    expand = c(0, 0),
    breaks = monthbreaks,
    labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

ggsave(file="bites/figs/CalendarYear.jpeg", width=10, height=2)
