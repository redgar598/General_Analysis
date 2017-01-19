###Explotory data analysis

library(rafalib)
data(father.son,package="UsingR") ##available from CRAN
x <- father.son$fheight

## see is heights fall along the normal distribution
ps <- ( seq(0,99) + 0.5 )/100 
qs <- quantile(x, ps)
normalqs <- qnorm(ps, mean(x), popsd(x))
plot(normalqs,qs,xlab="Normal percentiles",ylab="Height percentiles")
abline(0,1) ##identity line

library(ggplot2)
plt<-data.frame(normalqs=normalqs,qs=qs )
ggplot(plt, aes(normalqs,qs))+geom_point(color="cornflowerblue")+
  geom_abline(intercept = 0, slope = 1, color="grey40")+theme_bw()+
  xlab("Normal percentiles")+ylab("Height percentiles")



## non-normal data qq plot
# this data pulled from a t distribution which has fatter tails

dfs <- c(3,6,12,30)
mypar(2,2)
for(df in dfs){
  x <- rt(1000,df)
  qqnorm(x,xlab="t quantiles",main=paste0("d.f=",df),ylim=c(-6,6))
  qqline(x)
}

plt2<-lapply(dfs, function(df) {x<-rt(1000, df)
data.frame(t_quantiles=quantile(x, ps), Sample_quantiles=qnorm(ps, mean(x), popsd(x)),df=df)})
plt2<-do.call(rbind, plt2)

ggplot(plt2, aes(Sample_quantiles,t_quantiles))+geom_point(color="cornflowerblue")+
  geom_abline(intercept = 0, slope = 1, color="grey40")+theme_bw()+
  xlab("t quantiles")+ylab("Sample Quantiles")+facet_wrap(~df)



########### 
## BOX PLOTS
###########  
library(UsingR)
mypar(1,2)
hist(exec.pay) ##in UsingR package
qqnorm(exec.pay)
qqline(exec.pay)

boxplot(exec.pay, ylab="10,000s of dollars", ylim=c(0,400))

ggplot(as.data.frame(exec.pay), aes(1, exec.pay))+geom_violin(fill="grey", color="white")+geom_boxplot(width=0.2)+theme_bw()+ylim(0,400)


########### 
## Scatter plots and correlation
###########  

data("father.son")
x=father.son$fheight
y=father.son$sheight
plot(x,y,xlab="Father's height in inches",ylab="Son's height in inches",main=paste("correlation =",signif(cor(x,y),2)))

ggplot(father.son, aes(fheight, sheight))+geom_point(color="cornflowerblue")+theme_bw()+
  xlab("Father's height in inches")+ylab("Son's height in inches")

########### 
## Stratification
###########  
groups <- split(y,round(x))
boxplot(groups)
print(mean(y[ round(x) == 72]))

library(reshape2)
ggplot(melt(groups), aes(L1, value, fill=L1))+geom_boxplot()+theme_bw()

########### 
## Bi-variate Normal distribution
###########    
groups <- split(y,round(x))
mypar(2,2)
for(i in c(5,8,11,14)){
  qqnorm(groups[[i]],main=paste0("X=",names(groups)[i]," strata"),
         ylim=range(y),xlim=c(-2.5,2.5))
  qqline(groups[[i]])
}


library(psych)
scatter.hist(x=father.son$fheight, y=father.son$sheight, density=TRUE, ellipse=TRUE)

ggplot(melt(groups)[which(melt(groups)$L1%in%c(63,66,69,72)),], aes(L1, value))+
  geom_violin(fill="grey", color="white")+geom_boxplot(width=0.1)+theme_bw()+
  xlab("Father Height Strata")+ylab("Son Height")

#center/standardize
x=( x-mean(x) )/sd(x) #father.son$fheight
y=( y-mean(y) )/sd(y) #father.son$sheight
means=tapply(y, round(x*4)/4, mean)#rounds the son heights to the nearest 0.25 then startifies farther height by those intervals and takes mean
fatherheights=as.numeric(names(means))
mypar(1,1)
plot(fatherheights, means, ylab="average of strata of son heights", ylim=range(fatherheights))
abline(0, cor(x,y))

########### 
## Plots to avoid
########### 
library("downloader")
filename <- "fig1.RData"
url <- "https://github.com/kbroman/Talk_Graphs/raw/master/R/fig1.RData"
if (!file.exists(filename)) download(url,filename)
load(filename)

library(rafalib)
mypar()
dat <- list(Treatment=x,Control=y)
boxplot(dat,xlab="Group",ylab="Response",cex=0)
stripchart(dat,vertical=TRUE,method="jitter",pch=16,add=TRUE,col=1)

ggplot(melt(dat), aes(L1, value))+geom_boxplot(fill="lightgrey", width=0.5)+geom_point(color="cornflowerblue",size=3,position=position_jitter(w=0.25))+theme_bw()+xlab("Group")+ylab("Response")


url <- "https://github.com/kbroman/Talk_Graphs/raw/master/R/fig4.RData"
filename <- "fig4.RData"
if (!file.exists(filename)) download(url, filename)
load(filename)
mypar(1,2)
plot(x,y,lwd=2,type="n")
fit <- lm(y~x)
abline(fit$coef,lwd=2)
b <- round(fit$coef,4)
text(78, 200, paste("y =", b[1], "+", b[2], "x"), adj=c(0,0.5))
rho <- round(cor(x,y),4)
text(78, 187,expression(paste(rho," = 0.8567")),adj=c(0,0.5))
plot(x,y,lwd=2)
fit <- lm(y~x)
abline(fit$coef,lwd=2)

