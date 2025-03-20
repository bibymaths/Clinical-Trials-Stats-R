## ----mean medians, echo=FALSE-----------------------------------------------------------
library(pander)  
set.seed(1)
control = c(163.1, 174.2, 150.2, 183.4, 184.2, 166.6, 176.9, 174.7, 151.2, 174.9)
treat   = c(169.8, 164.9, 158.8, 171.1, 148.8, 164, 187.5, 161.6, 179.5, 153.0)  
mean1 <- mean(control) 
med1 <- median(control)
mean2 <- mean(treat) 
med2 <- median(control) 
d1 <- as.data.frame(cbind(rbind(mean1, mean2), rbind(med1, med2))) 
rownames(d1) <- c("Control", "Treatment") 
colnames(d1) <- c("Mean","Median") 
pander(d1, caption = "Mean and the Median for body weight from rats per group")



## ----sd, echo=FALSE---------------------------------------------------------------------
set.seed(1)
s1 <- sd(control)
s2 <- sd(treat) 
d2 <- as.data.frame(rbind(s1,s2)) 
rownames(d2) <- c("Control", "Treatment") 
colnames(d2) <- "Standard Deviation" 
pander(d2, caption = "Standard Deviation for body weight from rats per group")


## ----rat data, echo=FALSE---------------------------------------------------------------
d3 <- as.data.frame(cbind(control,treat)) 
colnames(d3) <- c("Control(X1)", "Treatment(X2)") 
pander(d3,caption = "Body weight from rats in a control (X1) and an active (X2) treatment group")


## ----data, echo=FALSE, fig.align='center', fig.cap="Boxplot for treatment and control group", message=FALSE, warning=FALSE, out.width='50%'----
library(rankFD)  
library(brunnermunzel)
d4 = data.frame(Weight=c(treat,control),  
                Treatment=factor(c(rep(1,10),rep(2,10))))  
boxplot(Weight ~ Treatment, data = d4, col = "white",  
        names = c("Treatment","Control"),  
        xlab = "Group", ylab = "Body Weight") 
stripchart(Weight ~ Treatment, data = d4, method = "jitter",  
           pch = 19, col = 2:4, vertical = TRUE, add = TRUE)


## ----rel, echo=TRUE, message=FALSE, warning=FALSE, results='markup'---------------------
set.seed(1)
test1 <- rank.two.samples(Weight ~ Treatment, data = d4 ,  
                          method = "t.app", shift.int = FALSE)
print(test1)


## ----brunmunz, echo=TRUE, message=FALSE, warning=FALSE, results='markup'----------------
brunnermunzel.test(treat, control)


## ----confplot, echo=TRUE, fig.align='center', fig.cap="Confidence interval for relative effects on body weight data from rats", message=FALSE, warning=FALSE, out.width='50%',tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(test1)


## ----simulations, echo = TRUE-----------------------------------------------------------
set.seed(1)  
## Two standard normal samples
a1 = rnorm(8, mean = 0, sd = 1) 
a2 = rnorm(4, mean = 0, sd = 1)   
## Two exponential samples
b1 = rexp(8, rate = 1) 
b2 = rexp(4, rate = 1)  
## A standard normal and an N(0,9) sample
c1 = rnorm(8, mean = 0, sd = 1) 
c2 = rnorm(4, mean = 0, sd = 9) 
nsim = 1000 
 
outcome1 <- c(a1,a2) 
outcome2 <- c(b1,b2)  
outcome3 <- c(c1,c2) 
group <- c(0,0,0,0,0,0,0,0,1,1,1,1) 

data1 <- data.frame(outcome1 ,group) 
data2 <- data.frame(outcome2 ,group)
data3 <- data.frame(outcome3 ,group) 
 
rank1 <- rank(outcome1,ties.method="average") 
rank2 <- rank(outcome2,ties.method="average") 
rank3 <- rank(outcome3,ties.method="average") 
 
r2wper1 <- rep(NA, nsim)  
r2wper2 <- rep(NA, nsim)  
r2wper3 <- rep(NA, nsim) 
  
## Save the value of the test statistic  
## (i.e. rank sum in the second sample)  
## for each simulation run  

for (i in 1:nsim){
r2wper1[i] <- sum(sample(rank1,size=12,replace=FALSE)[9:12]) 
r2wper2[i] <- sum(sample(rank2,size=12,replace=FALSE)[9:12]) 
r2wper3[i] <- sum(sample(rank3,size=12,replace=FALSE)[9:12]) 
} 
 
w1 = 31 
w2 = 19 
w3 = 36 
  
p_right1 <- mean(r2wper1 >= w1) 
p_left1 <- mean(r2wper1 <= w1) 
p_twosided1 <- 2*min(p_left1,p_right1) 
  
p_right2 <- mean(r2wper2 >= w2) 
p_left2 <- mean(r2wper2 <= w2) 
p_twosided2 <- 2*min(p_left2,p_right2) 
 
p_right3 <- mean(r2wper3 >= w3) 
p_left3 <- mean(r2wper3 <= w3) 
p_twosided3 <- 2*min(p_left3,p_right3) 
 
res1 <- table(r2wper1)/nsim 
res2 <- table(r2wper2)/nsim 
res3 <- table(r2wper3)/nsim   


## ----simplot1, echo =FALSE, fig.align='center', fig.cap="Approximation of the distribution of the exact WMW test statistic under H0 for two standard normal samples", message=FALSE, warning=FALSE, out.width='50%',tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(res1,las=1,ylab="Permutation pmf",xlab="Rank sum",main="
Approximation of the distribution of the exact WMW test statistic under H0\nTwo standard normal samples") 


## ----simplot2, echo =FALSE, fig.align='center', fig.cap="Approximation of the distribution of the exact WMW test statistic under H0 for two exponential samples", message=FALSE, warning=FALSE, out.width='50%',tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(res2,las=1,ylab="Permutation pmf",xlab="Rank sum",main="
Approximation of the distribution of the exact WMW test statistic under H0\nTwo exponential samples") 


## ----simplot3, echo =FALSE, fig.align='center', fig.cap="Approximation of the distribution of the exact WMW test statistic under H0 for a standard normal and a N(0,9) sample", message=FALSE, warning=FALSE, out.width='50%',tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(res3,las=1,ylab="Permutation pmf",xlab="Rank sum",main="
Approximation of the distribution of the exact WMW test statistic under H0\nstandard normal and N(0,9) sample")


## ----basicrel, echo = TRUE--------------------------------------------------------------
p1 <- 1/12*(mean(rank1[8+ 1:4]) - mean(rank1[1:8])) + 0.5  
p2 <- 1/12*(mean(rank2[8+ 1:4]) - mean(rank2[1:8])) + 0.5   
p3 <- 1/12*(mean(rank3[8+ 1:4]) - mean(rank3[1:8])) + 0.5 


## ----releffecttable, echo=FALSE, results='asis'-----------------------------------------
re <- as.data.frame(cbind(p1, p2, p3)) 
colnames(re) <- c("p(a)", "p(b)" , "p(c)") 
rownames(re) <- "Relative Effect" 
pander(re, caption = "Relative effect for sample scenario: a - two standard normal, b - two exponential, c - standard normal and N(0,9)")


## ----ci_forrelativeeffect, echo = TRUE--------------------------------------------------
set.seed(1) 
N = 12 
n1 = 8 
n2 = 4 
R1 = rank(a1) 
R2 = rank(a2)
sigma1 = sum((rank1[1:n1]-R1-mean(rank1[1:n1])+(n1+1)/2)^2) / ((N-n1)^2*(n1-1))
sigma2 = sum((rank1[n1+1:n2]-R2-mean(rank1[n1+1:n2])+(n2+1)/2)^2) / ((N-n2)^2*(n2-1)) 
TN = (p1 - 1/2) / sqrt(sigma1/n1 + sigma2/n2) 
df = (sigma1/n1 + sigma2/n2)^2 / ((sigma1^2/(n1^2*(n1-1))) + (sigma2^2/(n2^2*(n2-1))))
Lower = p1 - qt(0.975, df = df) * sqrt(sigma1/n1 + sigma2/n2)
Upper = p1 + qt(0.975, df = df) * sqrt(sigma1/n1 + sigma2/n2) 
c(Lower, Upper)


## ----newtreat, echo=TRUE----------------------------------------------------------------
treat1  = c(169.8, 164.9, 158.8, 171.1, 148.8, 164, 187.5, 161.6, 179.5, 153.0)
treat2  = c(228.5, 203.6, 227.3, 230.7, 193.9, 228.3, 229, 213, 225.6, 227.5, 186.3)  
d5 = data.frame(Weight=c(treat1, treat2, control),  
                Treatment=factor(c(rep(1,10),rep(2,11),rep(3,10))))  


## ----treat2, echo=FALSE-----------------------------------------------------------------
pander(summary(treat2), caption = "Five number summary for the new treatment")


## ----rat boxplot, echo=FALSE, fig.align='center', fig.cap="Boxplot for two treatments and control groups", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
boxplot(Weight ~ Treatment, data = d5, xlab = "Group",  
        names = c("Treatment 1", "Treatment 2", "Control"),  
        ylab = "Body Weight", col="white") 
stripchart(Weight ~ Treatment, data = d5, method = "jitter",  
           pch = 19, col = 2:4, vertical = TRUE, add = TRUE)


## ----pairwise rel1, echo=TRUE, message=FALSE, results='markup'--------------------------
library(nparcomp)  
set.seed(1)
p13 <- npar.t.test(Weight~Treatment,data=subset(d5, Treatment != "2"),  
            method="t.app", rounds=5, info = FALSE) 
summary(p13)


## ----pairwise rel2, echo=TRUE, message=FALSE, results='markup'--------------------------
set.seed(1)
p23 <- npar.t.test(Weight~Treatment,data=subset(d5, Treatment != "1"),  
            method="t.app", rounds=5, info = FALSE) 
summary(p23)


## ----pairwise result, echo=FALSE, results='asis'----------------------------------------
res.pair <- as.data.frame(rbind(p13$Analysis, p23$Analysis))   
rownames(res.pair) <- NULL
pander(res.pair, caption = "Pairwise relative effects")


## ----pairwise plot1, echo=FALSE, fig.align='center', fig.cap="CI for relative effect: Treatment 1 and Control", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(p13) 


## ----pairwise plot2, echo=FALSE, fig.align='center', fig.cap="CI for relative effect: Treatment 2 and Control", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(p23) 


## ----pairwise rank, echo=TRUE-----------------------------------------------------------
set.seed(1)
pairank <- nparcomp(Weight~Treatment,data=d5,type="Dunnett", info = FALSE, 
                    control="1",asy.method="mult.t", rounds=5) 
summary(pairank)


## ----pairwise plot3, echo=FALSE, fig.align='center', fig.cap="Simultaneous CI for relative effect: p(1,3) - Treatment 1 and Control, p(1,2) - Treatment 2 and Control", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(pairank) 


## ----pairwise rankresult, echo=FALSE, results='asis'------------------------------------
resrank.pair <- as.data.frame(pairank$Analysis)   
rownames(resrank.pair) <- NULL
pander(resrank.pair, caption = "Rank based relative effects")


## ----kruskal, echo = TRUE, results='markup'---------------------------------------------
set.seed(1) 
kw <- rankFD(Weight~Treatment, data=d5, effect="unweighted", hypothesis="H0F")  
print(kw)


## ----kruskal plot ,echo=FALSE, fig.align='center', fig.cap="CI for relative effects: 1 - Treatment 1, 2 - Treatment 2, 3- Control", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(kw)


## ----protein conc, echo = TRUE----------------------------------------------------------
## control group (A)  
## three dosage groups (B - D)
A = c(36.3, 11.1, 51.0, 9.8, 96.1, 245.8, 76.7, 88.4, 22.7, 123.7) 
B = c(74.9, 102.9, 170.6, 128.5, 149.2, 4.3, 125.1, 54.3, 161.0, 25.3)
C = c(20.9, 178.4, 50.9, 103.2, 61.5, 44.0, 60.0, 37.1, 39.9, 30.1)
D = c(43.4, 36.5, 80.1, 169.2, 31.3, 12.6, 23.8, 23.6, 17.6, 347.3) 
pc = data.frame(Group = c(rep("A", 10), rep("B", 10), rep("C", 10),  
                           rep("D", 10)), Concentration = c(A,B,C,D)) 


## ----prot boxplot, echo=FALSE, fig.align='center', fig.cap="Boxplot for blood protein concentrations in a control group(A) and in three dosage groups(B - D)", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff = 40)----
boxplot(Concentration ~ Group, data = pc, col = "white") 
stripchart(Concentration ~ Group, data = pc, method = "jitter",  
           pch = 19, col = 2:4, vertical = TRUE, add = TRUE)


## ----prot densityplot, echo=FALSE, fig.align='center', fig.cap="Distributional density for blood protein concentrations in a control group (A) and in three dosage groups (B - D)", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
library(ggplot2)
ggplot(pc, aes(x=Concentration, fill=Group)) + geom_density(alpha=.3) + theme_classic()


## ----prot histogram, echo=FALSE, fig.align = 'center', fig.cap="Histogram for blood protein concentrations in a control group (A) and in three dosage groups (B - D)", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
ggplot(pc, aes(x=Concentration, fill=Group)) + geom_histogram(bins = 50) + theme_classic() + facet_grid(Group ~ .)


## ----contrast, results='asis', echo=FALSE-----------------------------------------------
library(multcomp) 
set.seed(1)
n<-c(10,10,10,10) 
pander(contrMat(n, type = "Tukey"), caption = "Multiple Comparisons of Means: Tukey Contrasts")


## ----mctp1, echo =TRUE, results='markup'------------------------------------------------
library(nparcomp) 
set.seed(1) 
mctp1 <- mctp(Concentration~Group ,data=pc ,type="Tukey",asy.method="mult.t",  
                  effect = "weighted", rounds=5, info = FALSE)  
mctp2 <- mctp(Concentration~Group ,data=pc ,type="Tukey",asy.method="mult.t",  
                  effect = "unweighted", rounds=5, info = FALSE)   
mctp3 <- nparcomp(Concentration ~ Group ,data=pc, info = FALSE,
                  type="Tukey",asy.method="mult.t",rounds=5)
summary(mctp1) # Weighted Global pseudo ranks


## ----mctp2, echo =TRUE, results='markup'------------------------------------------------
summary(mctp2) # Unweighted Global pseudo ranks


## ----np releffect1, echo=FALSE, results='asis'------------------------------------------
library(knitr)
knitr::kable(mctp1$Analysis,  
             caption = "Non-parametric weighted relative effects: blood protein concentration using Global pseudo ranks") 


## ----np releffect2, echo=FALSE, results='asis'------------------------------------------
knitr::kable(mctp2$Analysis,  
             caption = "Non-parametric un-weighted relative effects: blood protein concentration using Global pseudo ranks") 


## ----np releffectp1, echo=FALSE, results='asis'-----------------------------------------
knitr::kable(mctp1$Overall,  
             caption = "p-value for non-parametric weighted relative effect p(B,C) using Global pseudo ranks") 


## ----np releffectp2, echo=FALSE, results='asis'-----------------------------------------
knitr::kable(mctp2$Overall,  
             caption = "p-value for non-parametric un-weighted relative effect p(B,C) using Global pseudo ranks") 


## ----mctp ci plot1 ,echo=FALSE, fig.align='center', fig.cap="Simultaneous CI for non-parametric weighted relative effects using Global pseudo ranks", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(mctp1)


## ----mctp ci plot2 ,echo=FALSE, fig.align='center', fig.cap="Simultaneous CI for non-parametric un-weighted relative effects using Global pseudo ranks", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(mctp2)


## ----mctp3, echo =TRUE, results='markup'------------------------------------------------
summary(mctp3) # pairwise rankings


## ----np releffect3, echo=FALSE, results='asis'------------------------------------------
knitr::kable(mctp3$Analysis,  
             caption = "Non-parametric relative effects: blood protein concentration using pairwise ranks") 


## ----np releffectp3, echo=FALSE, results='asis'-----------------------------------------
knitr::kable(mctp3$Overall,  
             caption = "p-value for non-parametric relative effect p(B,C) using pairwise ranks") 


## ----mctp ci plot3 ,echo=FALSE, fig.align='center', fig.cap="Simultaneous CI for non-parametric relative effects using pairwise ranks", message=FALSE, warning=FALSE, out.width='50%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
plot(mctp2)


## ----system, echo=FALSE-----------------------------------------------------------------
library(pander) 
pander(sessionInfo(), compact = FALSE)

