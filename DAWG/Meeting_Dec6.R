##### Data Analysis Working Group Meeting 2
setwd("~/Documents/Side/DAWG")

### Download data
library(dplyr)
library(rafalib)
library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- "mice_pheno.csv"
download(url,destfile=filename)

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if(!file.exists("femaleMiceWeights.csv")) download(url,destfile=filename)





###################
## Confidence intervals and t-tests
###################

dat <- read.csv("mice_pheno.csv")
chowPopulation <- dat[dat$Sex=="F" & dat$Diet=="chow",3]
mu_chow <- mean(chowPopulation)
print(mu_chow)

## can't get whole population so let's estimate!
N <- 30
chow <- sample(chowPopulation,N)
print(mean(chow))# pretty damn close
se <- sd(chow)/sqrt(N)
print(se)

## simulation to show CI (one example)
Q <- qnorm(1- 0.05/2)
interval <- c(mean(chow)-Q*se, mean(chow)+Q*se )
interval
interval[1] < mu_chow & interval[2] > mu_chow

## simulation to show CI (many examples to approach 95%)

B <- 250
mypar()
plot(mean(chowPopulation)+c(-7,7),c(1,1),type="n",
     xlab="weight",ylab="interval",ylim=c(1,B))
abline(v=mean(chowPopulation))
for (i in 1:B) {
  chow <- sample(chowPopulation,N)
  se <- sd(chow)/sqrt(N)
  interval <- c(mean(chow)-Q*se, mean(chow)+Q*se)
  covered <-
    mean(chowPopulation) <= interval[2] & mean(chowPopulation) >= interval[1]
  color <- ifelse(covered,1,2)
  lines(interval, c(i,i),col=color)
}


# Now with less samples!!
mypar()
plot(mean(chowPopulation)+c(-7,7),c(1,1),type="n",
     xlab="weight",ylab="interval",ylim=c(1,B))
abline(v=mean(chowPopulation))
Q <- qnorm(1- 0.05/2)
N <- 5
for (i in 1:B) {
  chow <- sample(chowPopulation,N)
  se <- sd(chow)/sqrt(N)
  interval <- c(mean(chow)-Q*se, mean(chow)+Q*se)
  covered <- mean(chowPopulation) <= interval[2] & mean(chowPopulation) >= interval[1]
  color <- ifelse(covered,1,2)
  lines(interval, c(i,i),col=color)
}

### this is worsened becasue we assumed a normal distribution when a t distribtuion may be more approporaite
mypar()
plot(mean(chowPopulation) + c(-7,7), c(1,1), type="n",
     xlab="weight", ylab="interval", ylim=c(1,B))
abline(v=mean(chowPopulation))
##Q <- qnorm(1- 0.05/2) ##no longer normal so use:
Q <- qt(1- 0.05/2, df=4)
N <- 5
for (i in 1:B) {
  chow <- sample(chowPopulation, N)
  se <- sd(chow)/sqrt(N)
  interval <- c(mean(chow)-Q*se, mean(chow)+Q*se )
  covered <- mean(chowPopulation) <= interval[2] & mean(chowPopulation) >= interval[1]
  color <- ifelse(covered,1,2)
  lines(interval, c(i,i),col=color)
}


####### EXCERCISE
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)

# birth weight smoking and non-smoking
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

# true population differences
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

#Q1
    set.seed(1)
    dat.ns<-bwt.nonsmoke[sample(1:length(bwt.nonsmoke), 25)]
    set.seed(1)
    dat.s<-bwt.smoke[sample(1:length(bwt.smoke), 25)]
    
    # standard error
    se <- sqrt(
    var(dat.s)/length(dat.s) +
      var(dat.ns)/length(dat.ns)
    )
    
    #tstat <- diff/se
    tval<-(mean(dat.s) - mean(dat.ns))/se

#Q2
    righttail <- 1 - pnorm(abs(tval))
    lefttail <- pnorm(-abs(tval))
    pval <- lefttail + righttail
    print(pval)
    
#Q3
    2*(pnorm(-abs(tval)))   #D
    
#Q4
    diff<-mean(dat.s) - mean( dat.ns)
    #99% confidence interval
    Q <- qnorm(1- 0.01/2) # critical value/z multiple/area under normal bounded at 1% and 99th %
    Q*se # amount to add and subtract from the mean to get a 99% confidence interval
    
            
            interval <- c(mean(diff)-Q*se, mean(diff)+Q*se )
            interval
            #In summary, if a 99% confidence interval does not include 0, then the p-value must be smaller than 0.01.
            
#Q5
    # with t distirbution to calculate Q instead of qnorm
    Q<-qt(0.01/2,df=2*25-2)
    Q*se
    
#Q6 C
    
#Q7 C
    
#Q8
    set.seed(1)
    dat.ns.mini<-bwt.nonsmoke[sample(1:length(bwt.nonsmoke), 5)]
    set.seed(1)
    dat.s.mini<-bwt.smoke[sample(1:length(bwt.smoke), 5)]
    
    pval<-t.test(dat.s.mini, dat.ns.mini)$p.value #,conf.level=0.9
    pval #0.02038
    
#Q9 B (logic on the previous statement of type II, power and alpha to rule out D)
    
#Q10
    set.seed(1)

    calculatePvalue <- function(N) {
      dat.ns.mini<-bwt.nonsmoke[sample(1:length(bwt.nonsmoke), N)]
      dat.s.mini<-bwt.smoke[sample(1:length(bwt.smoke), N)]
      pval<-t.test(dat.s.mini, dat.ns.mini)$p.value
    }

     pvalues <- sapply(rep(5, 10000), calculatePvalue)
    (length(which(pvalues<0.05))/10000)*100
    #9.84% of the time

#Q11
     set.seed(1)
     Ns <- seq(30,120,by=30)

     sapply(Ns, function(x){
       pvalues <- sapply(rep(x, 10000), calculatePvalue)
       (length(which(pvalues<0.05))/10000)*100
       })
     #N=60
     
#Q12
     set.seed(1)
     Ns <- seq(30,120,by=30)
     
     sapply(Ns, function(x){
       pvalues <- sapply(rep(x, 10000), calculatePvalue)
       (length(which(pvalues<0.01))/10000)*100}
     )
     #N=90
     
     
     
###################
## Monte Carlo Simulation
###################
library(dplyr)
dat <- read.csv("mice_pheno.csv")
controlPopulation <- filter(dat,Sex == "F" & Diet == "chow") %>%
 select(Bodyweight) %>% unlist

ttestgenerator <- function(n) {
  #note that here we have a false "high fat" group where we actually
  #sample from the nonsmokers. this is because we are modeling the *null*
  cases <- sample(controlPopulation,n)
  controls <- sample(controlPopulation,n)
  tstat <- (mean(cases)-mean(controls)) /
    sqrt( var(cases)/n + var(controls)/n )
  return(tstat)
}
ttests <- replicate(1000, ttestgenerator(10))
hist(ttests)

## generate control data with the population parameters
controls<- rnorm(5000, mean=24, sd=3.5)

#test the difference from controls
ttestgenerator <- function(n, mean=24, sd=3.5) {
  cases <- rnorm(n,mean,sd)
  controls <- rnorm(n,mean,sd)
  tstat <- (mean(cases)-mean(controls)) /
    sqrt( var(cases)/n + var(controls)/n )
  return(tstat)
}

####### EXCERCISE
#Q1
    set.seed(1)
    t=sqrt(5)*(mean(rnorm(5)))/sd(rnorm(5))
    t
#Q2
    set.seed(1)
    B<-sapply(rep(5, 1000), function(n) sqrt(n)*(mean(rnorm(n)))/sd(rnorm(n)))
    (length(which(B>2))/1000)*100
    
#Q3
    1-pt(2,df=4)
    B=100
    ps = seq(1/(B+1), 1-1/(B+1),len=B)
    qt(ps,df=4)
    
    qqplot(qt(ps,df=4),B,xlim=c(-4,4),ylim=c(-4,4))
    abline(0,1)
    
    ## now for several values of n not just 5
    sapply(seq(5,25,5), function(n){
    B<-sapply(rep(5, 1000), function(n) sqrt(n)*(mean(rnorm(n)))/sd(rnorm(n)))
    qqplot(qt(ps,df=n-1),B,xlim=c(-10,10),ylim=c(-10,10))
    abline(0,1)
    })
    #B?
    
#Q4
    ttestgenerator <- function(n) {
      cases <- rnorm(n)
      controls <- rnorm(n)
      tstat <- (mean(cases)-mean(controls)) /
        sqrt( var(cases)/n + var(controls)/n )
      return(tstat)
    }
    
    sapply(seq(10,100,10),function(Ns){
    ttests <- replicate(1000, ttestgenerator(Ns))
    
    ps <- (seq(0,999)+0.5)/1000
    qqplot(qt(ps,df=2*Ns-2),ttests,xlim=c(-10,10),ylim=c(-10,10))
    abline(0,1)
    })
    #C?
    
#Q5
    X=rbinom(n=15,size=1,prob=0.5)
    tstat <- sqrt(15)*mean(X) / sd(X)
   
    qqplot(qt(ps,df=14),X,xlim=c(-10,10),ylim=c(-10,10))
    abline(0,1)
    # False?
    
#Q6
    X=rbinom(n=500,size=1,prob=0.5)
    
    qqplot(qt(ps,df=499),X,xlim=c(-10,10),ylim=c(-10,10))
    abline(0,1)
    # False?
    
#Q7
    X=median(rnorm(500))
    B<-sapply(rep(500, 1000), function(n) median(rnorm(n)))
   
    #A)
    qqplot(rnorm(500, mean=0, sd=1/sqrt(500)),B)
    abline(0,1)
    
    #B)
    qqplot(rnorm(500),B)
    abline(0,1)
    
    #C)
    qqplot(qt(ps,df=499),B)
    abline(0,1)
    
    qqplot(qt(ps,df=4),sapply(rep(5, 1000), function(n) median(rnorm(n))))
    abline(0,1)
    
    #D)
    qqplot(rnorm(500, mean=0, sd=0.25),B)
    abline(0,1)
    
    #A?
    
    
###################
## Permutation Tests
################### 
    ## examples
    ## permutations 450k-Analysis/interpretation/pvalue distribution permutation and plotting
    ## Monte carlo 450k-Analysis/interpretation/Simulations_Overlap.R
    
    
dat=read.csv("femaleMiceWeights.csv")
library(dplyr)
control <- filter(dat,Diet=="chow") %>% select(Bodyweight) %>% unlist
treatment <- filter(dat,Diet=="hf") %>% select(Bodyweight) %>% unlist
obsdiff <- mean(treatment)-mean(control)

## shuffle labels
N <- 12
avgdiff <- replicate(1000, {
  all <- sample(c(control,treatment))
  newcontrols <- all[1:N]
  newtreatments <- all[(N+1):(2*N)]
  return(mean(newtreatments) - mean(newcontrols))
})
hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)

# calculate a permutation pvalue
#We add a 1 to the numerator and denominator to account for misestimation of the p-value (for more details see Phipson and Smyth, Permutation P-values should never be zero
(sum(abs(avgdiff) > abs(obsdiff)) + 1) / (length(avgdiff) + 1)



## now shuuffle with smaller sample
N <- 5
control <- sample(control,N)
treatment <- sample(treatment,N)
obsdiff <- mean(treatment)- mean(control)

avgdiff <- replicate(1000, {
  all <- sample(c(control,treatment))
  newcontrols <- all[1:N]
  newtreatments <- all[(N+1):(2*N)]
  return(mean(newtreatments) - mean(newcontrols))
})
hist(avgdiff)
abline(v=obsdiff, col="red", lwd=2)



### EXCERCISE
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- filter(babies, smoke==0) %>% select(bwt) %>% unlist
bwt.smoke <- filter(babies, smoke==1) %>% select(bwt) %>% unlist

#Q1
    #observed difference
    N=10
    set.seed(1)
    nonsmokers <- sample(bwt.nonsmoke , N)
    smokers <- sample(bwt.smoke , N)
    obs <- mean(smokers) - mean(nonsmokers)
    
    # one permutation
    dat <- c(smokers,nonsmokers)
    shuffle <- sample( dat )
    smokersstar <- shuffle[1:N]
    nonsmokersstar <- shuffle[(N+1):(2*N)]
    mean(smokersstar)-mean(nonsmokersstar)
    
    # 1000 permutations
    set.seed(1)
    avgdiff <- replicate(1000, {
      dat <- c(smokers,nonsmokers)
      shuffle <- sample( dat )
      smokersstar <- shuffle[1:N]
      nonsmokersstar <- shuffle[(N+1):(2*N)]
      return(mean(smokersstar)-mean(nonsmokersstar))
    })
    hist(avgdiff)
    abline(v=obs, col="red", lwd=2)
    
    (sum(abs(avgdiff) > abs(obs)) + 1) / (length(avgdiff) + 1)
    
    ##  0.05294705
    
    
#Q2
    ## now median
    N=10
    set.seed(1)
    nonsmokers <- sample(bwt.nonsmoke , N)
    smokers <- sample(bwt.smoke , N)
    obs <- median(smokers) - median(nonsmokers)
    
    # 1000 permutations
    set.seed(1)
    avgdiff <- replicate(1000, {
      dat <- c(smokers,nonsmokers)
      shuffle <- sample( dat )
      smokersstar <- shuffle[1:N]
      nonsmokersstar <- shuffle[(N+1):(2*N)]
      return(median(smokersstar)-median(nonsmokersstar))
    })
    hist(avgdiff)
    abline(v=obs, col="red", lwd=2)
    
    (sum(abs(avgdiff) > abs(obs)) + 1) / (length(avgdiff) + 1)
    
    #0.01798202
    
    
    
###################
## Association Tests
###################
disease=factor(c(rep(0,180),rep(1,20),rep(0,40),rep(1,10)),
               labels=c("control","cases"))
genotype=factor(c(rep("AA/Aa",200),rep("aa",50)),
                levels=c("AA/Aa","aa"))
dat <- data.frame(disease, genotype)
dat <- dat[sample(nrow(dat)),] #shuffle them up
head(dat)
table(genotype)
table(disease)
tab <- table(genotype,disease)
tab

## odds ratio (explain well in text)
(tab[2,2]/tab[2,1]) / (tab[1,2]/tab[1,1])

## calculate a p value to go with the OR need to build expectation 
p=mean(disease=="cases")
p
expected <- rbind(c(1-p,p)*sum(genotype=="AA/Aa"),
                  c(1-p,p)*sum(genotype=="aa"))
dimnames(expected)<-dimnames(tab)
expected
# chisq tests this on the tab of observed (builds its own expected)
chisq.test(tab)$p.value


## increase the N changes the p but no the OR
tab<-tab*10
chisq.test(tab)$p.value

## since OR is no normal (ratio of ratio) it is converted to log OR using GLM in order to build confidence intervals
fit <- glm(disease~genotype,family="binomial",data=dat)
coeftab<- summary(fit)$coef
coeftab

## caluclate confidence interval from the estimate and se
ci <- coeftab[2,1] + c(-2,2)*coeftab[2,2] # CI includes 1 and is therefore fail to reject null
exp(ci)

## EXCERCISE

#Q1
    url <- "https://studio.edx.org/c4x/HarvardX/PH525.1x/asset/assoctest.csv"
    filename <- basename(url)
    download(url, destfile=filename)
    d <- read.csv("assoctest.csv", header=TRUE)
    
    tab <- table(d$allele,d$case)
    tab
    chisq.test(tab)
    #X-squared = 3.3437
    
#Q2
    fisher.test(tab)
    #0.05194
