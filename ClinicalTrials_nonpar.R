##---------------------------------------------------------------
##                  Programme: Bioinformatics                   -
##                Matriculation number: 5491552                 -
##                    Author: Abhinav Mishra                    -
##--------------------------------------------------------------- 
 
## ----CRM 1, echo=TRUE, fig.align='center', fig.cap="Posterior Distributions with gamma(1,1) as prior at the start of the trial", message=FALSE, warning=FALSE, out.width='80%',tidy=TRUE, tidy.opts=list(width.cutoff=40)----

library(ggplot2)
library(rstan)
library(ggmcmc)
library(ggpubr)
library(cowplot) 
library(knitr)
 
##----------------------------------------------------------------
##      Modifying the stan code for priori beta = gamma(1,1)     -
##                          beta = rate                          -
##                         alpha = shape                         -
##                       syntax restricted                       -
##                  Bet = model parameter beta                   - 
##              EDITED SNIPPET: transformed data {} block        -
##---------------------------------------------------------------- 

smodel = "data {
      int<lower=1> K; // total number of dose levels 
      int<lower=0> N[K]; // total number of observations per dose
      int<lower=0> tox[K]; // number of tox per dose
      vector[K] skeleton;  
      } 
      
      transformed data {
      real<lower=0> alpha; // const. unmodeled param
      real<lower=0> beta; // const. unmodeled param
      alpha = 1.0;
      beta = 1.0;
      }
      
      parameters {
      real Bet;
      }
      
      transformed parameters {  
      real<lower=0,upper=1> p[K]; 
      
      for(k in 1:K){
      p[k] = pow(skeleton[k],exp(Bet));
      }
      
      }
      
      model {
      tox ~ binomial(N,p);
      Bet ~ gamma(alpha, beta);
      }
      "
model.crm.gamma <- stan_model(model_code = smodel)
 
# unnecesaary, already defined in  
# the stan model but okay
alpha <- 1 
beta <- 1 
tox <- c(0,0,0,0,0,0)
N <- c(0,0,0,0,0,0)
K <- length(tox) 
 
## provided in question 
skeleton <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5) 
 
set.seed(1)  

samples <- sampling(model.crm.gamma, 
                    data = list("tox"=tox, "N"=N, "K"=K,  
                                "skeleton"=skeleton,  
                                "alpha"=alpha, "beta"=beta), warmup = 2000,  
                    iter = 10000, chains = 4, cores = 1, 
                    thin = 1, refresh = 0)  
 
vis_CRM <- function(samples=samples) { 
   
  data <- ggs(samples) 
 
##---------------------------------------------------------------
##          Plot the distributions of the probabilities         -
##        of toxicities for each of the six dose levels         -
##--------------------------------------------------------------- 

plot_beta <- ggplot(data = data[data$Parameter=="Bet",], aes(x=value)) + geom_density(fill="dodgerblue3") +  
  xlab("Value") + ggtitle("beta~gamma(1,1)") + theme(axis.title.x = element_text(size=11),  
                                                                axis.title.y = element_text(size=11)) + 
  coord_cartesian(xlim = c(-6,6)) + theme_classic() 

plot1 <- ggplot(data = data[data$Parameter=="p[1]",], aes(x=value)) +  
  geom_density(fill="lightblue") + xlab("ProbTox") + ggtitle("Dose 1") + 
  coord_cartesian(ylim = c(0,9),  
                  xlim = c(0,1)) +  
  theme(axis.title.x = element_text(size=11),  
        axis.title.y = element_text(size=11),  
        plot.title = element_text(vjust=-7)) + theme_classic() 

plot2 <- ggplot(data = data[data$Parameter=="p[2]",], aes(x=value)) +  
  geom_density(fill="lightblue") + xlab("ProbTox") + ggtitle("Dose 2") +  
  coord_cartesian(ylim = c(0,9),  
                  xlim = c(0,1)) +  
  theme(axis.title.x = element_text(size=11),  
  axis.title.y = element_text(size=11),  
  plot.title = element_text(vjust=-7)) + theme_classic() 

plot3 <- ggplot(data = data[data$Parameter=="p[3]",], aes(x=value)) +  
  geom_density(fill="lightblue") + xlab("ProbTox") + ggtitle("Dose 3") +  
  coord_cartesian(ylim = c(0,9),  
                  xlim = c(0,1)) +  
  theme(axis.title.x = element_text(size=11),  
  axis.title.y = element_text(size=11),  
  plot.title = element_text(vjust=-7)) + theme_classic() 

plot4 <- ggplot(data = data[data$Parameter=="p[4]",], aes(x=value)) +  
  geom_density(fill="lightblue") + xlab("ProbTox") + ggtitle("Dose 4") +  
  coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11),  
                                                        axis.title.y = element_text(size=11),  
                                                        plot.title = element_text(vjust=-7)) + theme_classic() 

plot5 <- ggplot(data = data[data$Parameter=="p[5]",], aes(x=value)) +  
  geom_density(fill="lightblue") + xlab("ProbTox") + ggtitle("Dose 5") +  
  coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11),  
                                                        axis.title.y = element_text(size=11),  
                                                        plot.title = element_text(vjust=-7)) + theme_classic() 

plot6 <- ggplot(data = data[data$Parameter=="p[6]",], aes(x=value)) +  
  geom_density(fill="lightblue") + xlab("ProbTox") + ggtitle("Dose 6") +  
  coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11),  
                                                        axis.title.y = element_text(size=11),  
                                                        plot.title = element_text(vjust=-7)) + theme_classic() 

plots <- ggarrange(plot1, plot2, plot3, plot4,  
                   plot5, plot6, ncol=2, nrow=3) 

from_beta_to_prior_p <- ggdraw() +  
  draw_plot(plot_beta, x=0, y=0.25, width=0.4, height=0.5) +   
  draw_plot(plots, x=0.4, y=0, width=0.6, height=1) 

return(from_beta_to_prior_p)  
 
}
 
vis_CRM(samples=samples)


## ----CRM 2, echo=TRUE, message=FALSE, warning=FALSE,  fig.cap= "Posterior Distributions with gamma(1,1) as prior, after first cohort", fig.align='center', out.width='80%', tidy=TRUE, tidy.opts=list(width.cutoff=40)----
 
library(pander)
 
N1 <- c(3,0,0,0,0,0)  
  
set.seed(1) 
 
## Rerun the draws after first cohort
samples1 <-
  sampling(
    model.crm.gamma,
    data = list(
      "tox" = tox,
      "N" = N1,
      "K" = K,
      "skeleton" = skeleton,
      "alpha" = alpha,
      "beta" = beta),
    warmup = 2000,
    iter = 10000,
    chains = 4,
    cores = 1,
    thin = 1,
    refresh = 0
  )   
 
## plot the posterior distributions for each dose
stan_plot(samples1,  
          pars = c("p[1]","p[2]","p[3]","p[4]","p[5]","p[6]"),  
          ci_level = 0.95,  
          show_density = TRUE) +  
  ggtitle("Posterior Distributions for each dose (After 1st cohort)")
 
ssbeta <- rstan::summary(samples1, "Bet")$summary 
row.names(ssbeta)<-("beta") 
 
## report summary statistics for beta
knitr::kable(ssbeta, digits = 4,  
             caption = "Summary statistics for beta (After 1st cohort)")  
  

ci_p <- summary(samples1,  
                pars =c("p[1]","p[2]","p[3]","p[4]","p[5]","p[6]"),  
                probs = c(0.05, 0.95))$summary  

ci_95p <- ci_p[, c("5%","95%")]    

## 95% credible intervals for each dose level
pander(ci_95p, caption = "95% credible intervals for each dose level")



## ----CRM 3, echo=TRUE, message=FALSE, warning=FALSE,  fig.cap= "Posterior Distributions with gamma(1,1) as prior at the end of the trial", fig.align='center', out.width='80%'----
 
N2 <- c(3, 9, 9, 3, 0, 0)   
tox1 <- c(0, 1, 5, 1, 0, 0)
  
set.seed(1) 
 
## Rerun the draws after 24 patients
samples2 <-
  sampling(
    model.crm.gamma,
    data = list(
      "tox" = tox1,
      "N" = N2,
      "K" = K,
      "skeleton" = skeleton,
      "alpha" = alpha,
      "beta" = beta),
    warmup = 2000,
    iter = 10000,
    chains = 4,
    cores = 1,
    thin = 1,
    refresh = 0
  )   
 
## plot the posterior distributions for each dose
stan_plot(samples2, pars = c("p[1]","p[2]","p[3]", 
                             "p[4]","p[5]","p[6]"),  
          ci_level = 0.95, show_density = TRUE) +  
  ggtitle("Posterior Distributions  
          for each dose (After 24 patients)")
 
ssbeta1 <- rstan::summary(samples2, "Bet")$summary 
row.names(ssbeta1)<-("beta") 
 
## report summary statistics for beta
knitr::kable(ssbeta1, digits = 4, caption =  
               "Summary statistics for beta (After 24 patients)")  
  

ci_p24 <- summary(samples2, pars =c("p[1]","p[2]","p[3]","p[4]","p[5]","p[6]"),  
                  probs = c(0.05, 0.95))$summary  

ci_95p24 <- ci_p[, c("5%","95%")]    

## 95% credible intervals for each dose level after 24 patients
pander(ci_95p24,  
       caption = "95% credible intervals for  
       each dose level after 24 patients")



## ----beta prior, echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE-----------------------------

##-----------------------------------------------------------------------------
##              Using the snippet from third exercise sheet Q1:               -
##               changing the skeleton with the one given here,               -
##  and the distance from 0.4 will be used for 40% with the third dose d[3]   -
##----------------------------------------------------------------------------- 
dk <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5) 
sd <- sqrt(1.5)
abs.dist <- 1
set.seed(1) 
while(abs.dist>0.01){
  beta <- rnorm(n=1e6, mean=0, sd=sd)
  p3 <- dk[3]^exp(beta) 
  quantile.95 <- as.numeric(quantile(p3, probs = 0.95)) 
  dist <- 0.4-quantile.95 
  abs.dist <- abs(dist)
  sd <- sd-0.01
}
sd
quantile.95


## ----CRMSpinback, echo=TRUE, message=FALSE, warning=FALSE,  fig.align='center', out.width='80%', fig.cap="Posterior Distributions with known prior", tidy=TRUE, tidy.opts=list(width.cutoff=40)----
 
prmodel = "data {
      int<lower=1> K; // total number of dose levels 
      int<lower=0> N[K]; // total number of observations (so far) per dose level
      int<lower=0> tox[K]; // number of tox per dose
      vector[K] skeleton;
      real mu;
      real<lower=0> sigma;
      }
      
      
      parameters {
      real beta;
      }
      
      transformed parameters {
      real<lower=0,upper=1> p[K];
      
      for(k in 1:K){
      p[k] = pow(skeleton[k],exp(beta));
      }
      
      }
      
      model {
      tox ~ binomial(N,p);
      beta ~ normal(mu, sigma);
      }
      "
model.crm.normal <- stan_model(model_code = prmodel) 
 
mu <- 0 
sigma <- 0.564 
 
tox <- c(0,0,0,0,0,0)
N <- c(0,0,0,0,0,0)
K <- length(tox) 
 
## provided in question 
skeleton <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.5) 
 
## setting the seed
set.seed(1)   
samples <-
  sampling(
    model.crm.normal,
    data = list(
      "tox" = tox,
      "N" = N,
      "K" = K,
      "skeleton" = skeleton,
      "mu" = mu,
      "sigma" = sigma
    ),
    warmup = 2000,
    iter = 10000,
    chains = 4,
    cores = 1,
    thin = 1,
    refresh = 0
  )

vis_CRM <- function(samples=samples) {data <- ggs(samples)
plot_beta <- ggplot(data = data[data$Parameter=="beta",], aes(x=value)) + geom_density(fill="dodgerblue3") + 
  xlab("Value") + ggtitle("beta") + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11)) + 
  coord_cartesian(xlim = c(-6,6)) + theme_classic()
plot1 <- ggplot(data = data[data$Parameter=="p[1]",], aes(x=value)) + geom_density(fill="lightblue") + 
  xlab("Prob. of tox.") + ggtitle("Dose 1") + coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11), plot.title = element_text(vjust=-7)) + theme_classic()
plot2 <- ggplot(data = data[data$Parameter=="p[2]",], aes(x=value)) + geom_density(fill="lightblue") + 
  xlab("Prob. of tox.") + ggtitle("Dose 2") + coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11), plot.title = element_text(vjust=-7)) + theme_classic()
plot3 <- ggplot(data = data[data$Parameter=="p[3]",], aes(x=value)) + geom_density(fill="lightblue") + 
  xlab("Prob. of tox.") + ggtitle("Dose 3") + coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11), plot.title = element_text(vjust=-7)) + theme_classic()
plot4 <- ggplot(data = data[data$Parameter=="p[4]",], aes(x=value)) + geom_density(fill="lightblue") + 
  xlab("Prob. of tox.") + ggtitle("Dose 4") + coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11), plot.title = element_text(vjust=-7)) + theme_classic()
plot5 <- ggplot(data = data[data$Parameter=="p[5]",], aes(x=value)) + geom_density(fill="lightblue") + 
  xlab("Prob. of tox.") + ggtitle("Dose 5") + coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11), plot.title = element_text(vjust=-7)) + theme_classic()
plot6 <- ggplot(data = data[data$Parameter=="p[6]",], aes(x=value)) + geom_density(fill="lightblue") + 
  xlab("Prob. of tox.") + ggtitle("Dose 6") + coord_cartesian(ylim = c(0,9), xlim = c(0,1)) + theme(axis.title.x = element_text(size=11), axis.title.y = element_text(size=11), plot.title = element_text(vjust=-7)) + theme_classic()
plots <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, ncol=2, nrow=3) 

from_beta_to_prior_p <- ggdraw() + draw_plot(plot_beta, x=0, y=0.25, width=0.4, height=0.5) + draw_plot(plots, x=0.4, y=0, width=0.6, height=1)
return(from_beta_to_prior_p) }
vis_CRM(samples=samples) 
 
prbeta <- rstan::summary(samples, "beta")$summary 
 
stan_plot(samples,ci_level = 0.95, show_density = TRUE) +  
  ggtitle("Posterior Distributions  
          for each dose (using prior approximation)") 
 
medians <- as.data.frame(as.matrix(summary(samples)$summary[ , "50%"]))  

row.names(medians) <- c("beta","Dose 1","Dose 2","Dose 3", 
                        "Dose 4","Dose 5","Dose 6", "Intercept") 
colnames(medians) <- "Prob(MTD)" 
 
## report summary statistics for beta
knitr::kable(prbeta, digits = 3, caption = "Summary statistics  
             for beta (using prior approximation)")   



## ----mtd, echo=TRUE-------------------------------------------------------------------------------------
knitr::kable(medians,digits = 2,caption="Probabilities of toxicity for each dose")


## ----metrics, echo=TRUE, message=FALSE, warning=FALSE, paged.print=FALSE--------------------------------
 
## Given dataset
n <- 230 
score <- c(5,4,3,2,1) 
pos <- c(50,30,15,3,2) 
neg <- c(6,10,15,44,55) 
brc <- as.data.frame(cbind(score,pos,neg)) 
colnames(brc) <- c("Score", "Present", "Absent")

##----------------------------------------------------------------
##                    2 x 2 contigency table                     -
##                          Cutoff = 2                           -
##----------------------------------------------------------------
Sum1 <- c(126,104,n) 
Absent <- c(31,99,130) 
Present <- c(95,5,100) 
 
cont <- as.data.frame(cbind(Present, Absent, Sum1)) 
colnames(cont) <- c("Disease Present", "Disease Absent","") 
row.names(cont) <- c("Test +ve","Test -ve","") 

##---------------------------------------------------------------
##            Specificity & Sensitivity (+95% CI)               - 
##            NPV & PPV (+95% CI)                               - 
##            Prevelance
##---------------------------------------------------------------
z <- 1.96
sens <- cont[1,1]/cont[3,1]
spec <- cont[2,2]/cont[3,2] 
ppv <- cont[1,1]/cont[1,3] 
npv <- cont[2,2]/cont[2,3] 
prev <- cont[3,1]/n
 
spec95 <- c(spec-z*sqrt((spec*(1-spec))/(cont[3,2])), 
            spec+z*sqrt((spec*(1-spec))/(cont[3,2]))) 
sens95 <- c(sens-z*sqrt((sens*(1-sens))/(cont[3,1])), 
            sens+z*sqrt((sens*(1-sens))/(cont[3,1])))  
npv95 <- c(npv-z*sqrt((npv*(1-npv))/(cont[2,3])),  
           npv+z*sqrt((npv*(1-npv))/(cont[2,3]))) 
ppv95 <- c(ppv-z*sqrt((ppv*(1-ppv))/(cont[1,3])),  
           ppv+z*sqrt((ppv*(1-ppv))/(cont[1,3])))  
  
b <- as.data.frame(rbind(sens, spec, prev, ppv, npv))  

row.names(b) <- c("Sensitivity", 
                  "Specificity",  
                  "Prevalence", 
                  "Positive Predicted Value", 
                  "Negative Predicted Value")  
colnames(b) <- "Score"  
 
cd <- as.data.frame(rbind(t(as.data.frame(sens95)), t(as.data.frame(spec95)), 
                          t(as.data.frame(ppv95)),t(as.data.frame(npv95)))) 
colnames(cd) <- c("Lower","Upper")
row.names(cd) <- c("Senstivity","Specificity", 
                  "Positive Predicted Value","Negative Predicted Value")

##----------------------------------------------------------------
##                    2 x 2 contigency table                     -
##                          Cutoff = 3                           -
##----------------------------------------------------------------
Sum2 <- c(96,134,n) 
Absent <- c(16,114,130) 
Present <- c(80,20,100) 
 
cont1 <- as.data.frame(cbind(Present, Absent, Sum2)) 
colnames(cont1) <- c("Disease Present", "Disease Absent", "") 
row.names(cont1)<-c("Test +ve","Test -ve","") 

sens1 <- cont1[1,1]/cont1[3,1]
spec1 <- cont1[2,2]/cont1[3,2] 
 
##----------------------------------------------------------------
##                    2 x 2 contigency table                     -
##                          Cutoff = 4                           -
##----------------------------------------------------------------
Sum3 <- c(56,174,n) 
Absent <- c(6,124,130) 
Present <- c(50,50,100) 
 
cont2 <- as.data.frame(cbind(Present, Absent, Sum3)) 
colnames(cont2) <- c("Disease Present", "Disease Absent", "") 
row.names(cont2)<-c("Test +ve","Test -ve","") 
 
sens2 <- cont2[1,1]/cont2[3,1]
spec2 <- cont2[2,2]/cont2[3,2]  

##----------------------------------------------------------------
##              ROC plot using values from 3 cutoffs             -
##---------------------------------------------------------------- 

senss <- c(sens, sens1, sens2) 
specs <- c(spec, spec1, spec2)
rocx <- as.data.frame(cbind(c(2,3,4),senss,specs)) 
colnames(rocx) <- c("Cutoff", "Sensitivity", "Specificity")



## ----metric_a, echo=FALSE-------------------------------------------------------------------------------
knitr::kable(cont, caption = "2x2 table (cutoff = 2)") 


## ----metric_bc, echo=FALSE------------------------------------------------------------------------------
knitr::kable(b, digits = 3, caption = "Metrics (cutoff = 2)")


## ----metric_bcd, echo=FALSE-----------------------------------------------------------------------------
knitr::kable(cd, digits = 3,  
             caption = "Metrics with 95% CIs (cutoff = 2)") 


## ----metric_e, echo=FALSE, fig.align='center', fig.cap= "ROC curve (for cutoffs=2,3,4)", out.width='60%', fig.show='hold'----
knitr::kable(rocx, caption = "Plotted values for ROC")   
plot((1-specs),senss ,xlab = "1-Specificity", 
     ylab = "Sensitivity",type = "s", main = "ROC")


## ----metric_f, echo=FALSE-------------------------------------------------------------------------------
knitr::kable(cont1, caption = "2x2 table (cutoff = 3)")   
knitr::kable(cont2, caption = "2x2 table (cutoff = 4)")  


## ----auc, echo=TRUE-------------------------------------------------------------------------------------
AUC <- round(sum(specs[1:3]*diff(c(0,1 - senss[1:3]))),2) 
AUC


## ----design, echo=TRUE, message=FALSE, warning=FALSE, results='asis'------------------------------------

library(rpact) 
set.seed(1) 

design90p <-
  getDesignGroupSequential(
    sided = 2,
    alpha = 0.05,
    beta = 0.1, 
    normalApproximation = TRUE,
    typeOfDesign = "P")  
 
design80p <-
  getDesignGroupSequential(
    sided = 2,
    alpha = 0.05,
    beta = 0.2,  
    normalApproximation = TRUE,
    typeOfDesign = "P") 
 
design90p50 <-
  getDesignGroupSequential(
    sided = 2,
    alpha = 0.05,
    beta = 0.1, 
    informationRates = c(0.5, 1), 
    normalApproximation = TRUE,
    typeOfDesign = "P") 

sample90p <- getSampleSizeMeans(
  design90p,
  alternative = 4,
  stDev = 2.5,
  allocationRatioPlanned = 1)  
  
sample80p <- getSampleSizeMeans(
  design80p,
  alternative = 4,
  stDev = 2.5,
  allocationRatioPlanned = 1)   
 
powerSD1 <- getPowerMeans(design90p,
  alternative = 4,  
  stDev = 2,  
  sided = 2,
  allocationRatioPlanned = 1,  
  maxNumberOfSubjects = 20,  
  alpha = 0.05) 
 
powerSD2 <- getPowerMeans(design90p,
  alternative = 4,  
  stDev = 3,  
  sided = 2,
  allocationRatioPlanned = 1,  
  maxNumberOfSubjects = 20,  
  alpha = 0.05)
 
sample90p50 <- getSampleSizeMeans(
  design90p50,
  alternative = 4,
  stDev = 2.5,
  allocationRatioPlanned = 1) 


## ----designa, echo=FALSE, results='asis'----------------------------------------------------------------
kable(summary(sample90p)) 


## ----designb, echo=FALSE, results='asis'----------------------------------------------------------------
kable(summary(sample80p)) 


## ----designc1, echo=FALSE, fig.cap="Boundaries: Power and Standard Deviation", fig.show='hold', message=FALSE, warning=FALSE, out.width='50%'----
plot(powerSD1, type = 1) 
plot(powerSD2, type = 1) 


## ----designc2, echo=FALSE, fig.cap="Boundaries Effect Scale: Power and Standard Deviation", fig.show='hold', message=FALSE, warning=FALSE, out.width='50%'----
plot(powerSD1, type = 2) 
plot(powerSD2, type = 2) 


## ----designc, echo=FALSE, results='asis'----------------------------------------------------------------
kable(summary(powerSD1))  
kable(summary(powerSD2)) 


## ----designd, echo=TRUE,results='asis'------------------------------------------------------------------
kable(summary(sample90p50))


## ----kaplanmeier, echo=TRUE, message=FALSE, warning=FALSE, fig.cap= "Kaplan-Meier curves with 95% CIs", tidy=TRUE, tidy.opts=list(width.cutoff=40)----
library(survival) 
library(survRM2) 
library(ggfortify)  


load("~/Documents/Freie/Methods for Clinical Trials/colon-cancer.Rdata") 
set.seed(1)
model_fit <- survfit(Surv(timeD/365, statusD) ~ rx , data = d)  

autoplot(model_fit) + 
 labs(x = "\n Survival Time (Years) ", y = "Survival Probabilities \n", 
      title = "Kaplan-Meier Plot" ) +  
 theme(plot.title = element_text(hjust = 0.5), 
 axis.title.x = element_text(face="bold", colour="#FF7A33", size = 12),
 axis.title.y = element_text(face="bold", colour="#FF7A33", size = 12),
 legend.title = element_text(face="bold", size = 10)) 



## ----median, echo=TRUE, message=FALSE, warning=FALSE----------------------------------------------------
pander(model_fit, caption = "Median survival times with 95% CIs")


## ----median7, echo=TRUE, message=FALSE, warning=FALSE---------------------------------------------------
pander(summary(model_fit, times = 7))


## ----logrank, echo=TRUE, message=FALSE, warning=FALSE---------------------------------------------------
lr <- survdiff(Surv(timeD/365, statusD) ~ rx , data = d) 
pander(lr) 


## ----cox, echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------
effsize <- coxph(Surv(timeD/365, statusD) ~ rx , data = d) 
pander(effsize)


## ----coxassump, echo=FALSE, fig.cap="Schoenefeld Residuals", fig.show='hold', message=FALSE, warning=FALSE, out.width='50%'----
cz <- cox.zph(coxph(Surv(timeD/365, statusD) ~ rx + age , data = d)) 
par(mar = c(4, 4, .1, .1))
plot(cz)


## ----coxassump1, echo=TRUE------------------------------------------------------------------------------
cz


## ----system, echo=FALSE---------------------------------------------------------------------------------
library(pander) 
pander(sessionInfo(), compact = FALSE)

