#all models  - no x transformation
setwd("C:/Users/gricardo/OneDrive - Australian Institute of Marine Science/R/2019 demo")
#data1 = read.table("https://pastebin.com/raw/v7Uy4iXs", header= TRUE,dec=",")  #hav4
data1 = read.table("https://pastebin.com/raw/zfrUha88", header= TRUE,dec=",")
#data1 = read.table("https://pastebin.com/raw/BaCAP3Sr", header= TRUE,dec=",")
data1

options(scipen = 999)  # turn off scientific notation

# Load packages (if not installed Tools->Install Packages)
# library(lattice)
# library(MASS)  
# library (rstan) #Bayesian models


##########################################
#Data organising
data1$raw.x <- as.numeric(as.character(data1$raw.x))
data1$suc <- as.numeric(as.character(data1$suc))
data1$tot <- as.numeric(as.character(data1$tot))
data1$fail = data1$tot  - data1$suc
data1$fail <- as.numeric(as.character(data1$tot))

data1$obs <- factor(formatC(1:nrow(data1), flag="0", width = 3))# unique tank ID for later on
min.x <- min(data1$raw.x)
max.x <- max(data1$raw.x)
lg.x <- 100
df.x <- data.frame(raw.x = seq(min.x, max.x, length = lg.x)) #setting up  new  data frame (df) defining log.x values to run 
vec.x =df.x[,1]
par(mfrow = c(1, 1)) #setting space for plots

########################################################################################
#1) Binomial GLM 
md1 <- glm(cbind(suc,(tot - suc)) ~ raw.x, family = binomial (link = logit),data = data1)
summary(md1)
sum(residuals(md1, type = "deviance")^2)/md1$df.residual # >1.5, overdispersed i.e greater variability in the data than would be expected from the statistical model  (Logan?)
sum(resid(md1, type = "pearson")^2) / (nrow(data1) - length(coef(md1))) #(Zuur)
#plot(md1)

GLMfun = function(i){
  mm <- model.matrix(~raw.x, df.x)  # build model matrix to sub the parameters into the equation
  eta <- mm %*% coef(i)  #sub in parameters to get mean on logit scale
  prediction  <- as.vector(exp(eta) / (1 + exp(eta)))  #back-transform mean from logit scale
  se    <- sqrt(diag(mm %*% vcov(i) %*% t(mm)))  #work out standard errors
  upper  <- exp(eta + 1.96 *se) /(1 + exp(eta  + 1.96 *se))  #work out upper 95% CI
  lower  <- exp(eta - 1.96 *se) /(1 + exp(eta  - 1.96 *se)) #work out lower 95% CI
  df = data.frame(vec.x, prediction, lower,upper)  #put everything in a data.frame
  return(df)
} #coef
df1=GLMfun(md1)

plot(data1$raw.x,(data1$suc / data1$tot),log = 'x', main="Binomial GLM") #first plot
lines(df1$vec.x, df1$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df1$vec.x, df1$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df1$vec.x, df1$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

########################################################################################
#2) Run a quasibinomial GLM
md2 <- glm(cbind(suc,(tot - suc)) ~ raw.x, family = quasibinomial (link = logit),data = data1)
summary(md2)
df2=GLMfun(md2)

plot(data1$raw.x,(data1$suc / data1$tot),log = 'x', main="Quasibinomial GLM") #second plot
lines(df2$vec.x, df2$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df2$vec.x, df2$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df2$vec.x, df2$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

#######################################################################################################
#3) Run a GLMM with a random slope effect (i.e tank / container etc)
library(lme4)
md3 <- glmer(cbind(suc,(tot - suc)) ~ raw.x+(1|obs),family = binomial (link = logit),data = data1) 
library(RVAideMemoire) #GLMM overdispersion test
overdisp.glmer(md3) #Overdispersion for GLMM

GLMfun2 = function(i){
  mm <- model.matrix(~raw.x, df.x)  # build model matrix 
  eta <- mm %*% fixef(i)
  prediction  <- as.vector(exp(eta) / (1 + exp(eta)))
  se    <- sqrt(diag(mm %*% vcov(i) %*% t(mm)))
  upper  <- as.vector(exp(eta + 1.96 *se) /(1 + exp(eta  + 1.96 *se)))
  lower  <- as.vector(exp(eta - 1.96 *se) /(1 + exp(eta  - 1.96 *se)))
  df = data.frame(vec.x, prediction, lower,upper)
  return(df)
} #fixef
df3 = GLMfun2(md3)

#plot
plot(data1$raw.x,(data1$suc / data1$tot),log = 'x',main="GLMM") #third plot
lines(df3$vec.x, df3$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df3$vec.x, df3$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur upper 95% CI
lines(df3$vec.x, df3$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur lower 95% CI


#########################################################
library(glmmTMB)
md3.beta <- glmmTMB(cbind(suc,(tot - suc)) ~ raw.x +(1|obs) , family='betabinomial', data= data1)
summary(md3.beta)


  mm <- model.matrix(~raw.x, df.x)  # build model matrix to sub the parameters into the equation
  eta <- mm %*% fixef(md3.beta)$cond  #sub in parameters to get mean on logit scale
  prediction  <- as.vector(exp(eta) / (1 + exp(eta)))  #back-transform mean from logit scale
  se    <-  sqrt(diag(mm %*% vcov(md3.beta)$cond %*% t(mm)))
  upper  <- exp(eta + 1.96 *se) /(1 + exp(eta  + 1.96 *se))  #work out upper 95% CI
  lower  <- exp(eta - 1.96 *se) /(1 + exp(eta  - 1.96 *se)) #work out lower 95% CI
  df = data.frame(vec.x, prediction, lower,upper)  #put everything in a data.frame
  

  plot(data1$raw.x,(data1$suc / data1$tot),log = 'x',main="betabin") #third plot
  lines(df$vec.x, df$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
  lines(df$vec.x, df$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur upper 95% CI
  lines(df$vec.x, df$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur lower 95% CI
  
  
  
  

############################################################################################
# 4) Non-linear regression (4-par logistic)
library(drc) #NLR
md4 <- drm((suc / tot)~ raw.x, fct = LL.4(fixed = c(NA, 0, NA, NA),names = c("Slope", "Lower", "Upper", "ED50")),data = data1) 

pred.m4 <- predict(md4, newdata = df.x, interval="confidence") #predict 95% CI
df4 = data.frame(vec.x, pred.m4)

plot(data1$raw.x,(data1$suc/ data1$tot),log="x", main="LL.4 Lower=0") #fouth plot
lines(df4$vec.x, df4$Prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df4$vec.x, df4$Lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df4$vec.x, df4$Upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

############################################################################################################
# 5) Non-linear regression (Weighted binomial 4-par logistic, (refered to as Generalised nonlinear models  - a logistic regression model with a parameter-dependant link function)
library(drc) #NLR
md5 <- drm((suc / tot)~ raw.x, fct = LL.4(fixed = c(NA, 0, NA, NA),names = c("Slope", "Lower", "Upper", "ED50")),  weights = tot, type = 'binomial',data = data1) 
#plot(fitted(md5), residuals(md5, typeRes = "standard"))

pred.m5 <- predict(md5, newdata = df.x, interval="confidence") #predict 95% CI
df5 = data.frame(vec.x, pred.m5)

plot(data1$raw.x,(data1$suc/ data1$tot),log="x", main="LL.4 Binomial Lower=0") #fouth plot
lines(df5$vec.x, df5$Prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df5$vec.x, df5$Lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df5$vec.x, df5$Upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

##################################################################################
# 6) Non-linear regression (Weighted binomial 4-par Weibull)
library(drc) #NLR
md6 <- drm((suc / tot)~ raw.x, fct = W1.4(fixed = c(NA, 0, NA, NA),names = c("Slope", "Lower", "Upper", "ED50")),data = data1) 

pred.m6 <- predict(md6, newdata = df.x, interval="confidence") #predict 95% CI
df6 = data.frame(vec.x, pred.m6)

plot(data1$raw.x,(data1$suc/ data1$tot),log="x", main="W1.4 Weibull Lower=0") #fouth plot
lines(df6$vec.x, df6$Prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df6$vec.x, df6$Lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df6$vec.x, df6$Upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

#anova(md4,md6) #diff. between models

# ########################################################################################
#7) Bayesian binomial GLM
#compile rtools code (avoid the 127 error!)
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

library (brms) #Bayesian models
md7 <- brm ( suc | trials(tot)~ raw.x, data=data1,
             family="binomial"(link="logit"), chains=3, iter=2000, warmup=500,
             prior=c(set_prior ("cauchy (0, 2.5)", class="b"))) #mean, SD (cauchy Gelman's weakly informative prior)
plot(md7)

GLMfun3 = function(i){
  mm <- model.matrix(~raw.x, df.x)  # build model matrix
  med.coef = c(fixef(i)[1, 1], fixef(i)[2, 1])
  eta <- mm %*% med.coef
  prediction  <- as.vector(exp(eta) / (1 + exp(eta)))
  se    <- sqrt(diag(mm %*% vcov(i) %*% t(mm)))
  upper  <- as.vector(exp(eta + 1.96 *se) /(1 + exp(eta  + 1.96 *se)))
  lower  <- as.vector(exp(eta - 1.96 *se) /(1 + exp(eta  - 1.96 *se)))
  df = data.frame(vec.x, prediction, lower,upper)
  return(df)
} #fixef
df7 = GLMfun3(md7)

plot(data1$raw.x,(data1$suc/data1$tot), log = 'x', main="Bayesian Binom GLM") #first plot
lines(df7$vec.x,df7$prediction, type = "l", lwd = 2, col=2, xaxt = "n", las = 1) #plot model mean line
lines(df7$vec.x,df7$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df7$vec.x,df7$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI
 
##############################################################################################
# #Bayesian binomial GLMM
library(brms)
md8 <- brm ( suc | trials(tot)~ raw.x +(1|obs), data=data1,
             family="binomial"(link="logit"), chains=3, iter=2000, warmup=500,
             prior=c(set_prior ("cauchy (0, 2.5)", class="b")))
plot(md8)

df8 = GLMfun3(md8)

plot(data1$raw.x,(data1$suc/data1$tot),log = 'x',main="Bayesian Binom GLMM") #first plot
lines(df8$vec.x,df8$prediction, type = "l", lwd = 2, col=2, xaxt = "n", las = 1) #plot model mean line
lines(df8$vec.x,df8$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df8$vec.x,df8$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

##################################################################################
#ANOVA NOEC  (don't do this!)
prop.t = asin(data1$suc/data1$tot)
md9 = aov(prop.t~factor(raw.x), data = data1)
summary(md9)
library(emmeans)
options(scipen = 999)  # turn off scientific notation
noec.ano = as.data.frame(contrast(lsmeans(md9, "raw.x") , "trt.vs.ctrl1", ref = "0.1"))
noec.ano
pos = max(which(noec.ano$ p.value >0.05))
noec.ano$contrast[pos]  #this code is dodgy, check contrast!
####################################################
#GLMM NOEC
raw.x.cat <- as.factor(as.character(data1$raw.x))
md10 <- glmer(cbind(suc,(tot - suc)) ~ factor(raw.x)+(1|obs),family = binomial (link = logit),data = data1) 
summary(md10)
library(emmeans)
noec.glmm = as.data.frame(contrast(lsmeans(md10, "raw.x") , "trt.vs.ctrl1", ref = "0.1"))
noec.glmm
pos = max(which(noec.glmm$ p.value >0.05))
noec.glmm$contrast[pos]


#########################################################
#Bayesian NEC

sim.dat.data1 <- as.list(data1) #Bugs loves lists
data1.data <- list(raw.x      = data1$raw.x,      #Response
                  suc = data1$suc, #Covariate
                  tot       = data1$tot,
                  N       = nrow(data1))  #Sample size

#Use this to visualise priors (you will need to checnge this in the model)
SD = 100
Reciprical  = 1/SD #recip of used in jags
1/0.001
plot(seq(-250,250,.5),dnorm(seq(-250,250,.5),0.9,(SD)),type="l") #plot distribution of the normal. Note: jags distributions are often the reciprical of the R distribution, so an SD of 1000 would be 0.0001
plot(seq(0,500,.5),dgamma(seq(0,500,.5),0.0001,0.0001),type="l")
plot(seq(0,1,.1),dnorm(seq(0,1,.1),0.9,1000),type="l")  #truncated

sink("data1_model.txt")
cat("
    model
    {
    # specify model priors
    top~dnorm(0.9,0.0001)T(0,1)
    beta~dgamma(0.0001,0.0001)
    NEC~dnorm(30, 0.0001) T(0,)   #mean set to 30, contrained to >0
    
    # likelihood
    for (i in 1:N)
    {
    # theta is true proportion at each dose - given by equation (2) in text
    theta[i]<-top*exp(-beta*(raw.x[i]-NEC)*step((raw.x[i]-NEC)))
    # response is suc (number surviving) - assumed to be binomial
    suc[i]~dbin(theta[i],tot[i])
    }
    }
    ", fill=TRUE)
sink()  #Made model in working directory

params <- c("top", "beta", "NEC")
inits <- function(){list(top = 0.9, beta = 0.2, NEC = 30)}  #you will need to change this

library(R2jags)
J1.data1 <- jags(data       = data1.data,
                inits      = inits,
                parameters = params,
                model      = "data1_model.txt",
                n.thin     = 10,
                n.chains   = 3,
                n.burnin   = 5000,
                n.iter     = 50000)

J2.data1  <- update(J1.data1, n.iter = 10000, n.thin = 10)  
out.data1 <- J2.data1$BUGSoutput
Beta.mcmc2 = out.data1$sims.list$NEC
source(file = "MCMCSupportHighstatV4.R")
source(file = "HighstatLibV10.R")
library(MASS) 
library(rjags)
library(coda)
library(lattice)
MyBUGSChains(out.data1,c("top", "beta", "NEC"))
print(out.data1, digits = 3)  #summary

#plotting
plot(data1$raw.x,(data1$suc/data1$tot),log = 'x',main="NEC-bayes") #first plot
NEC.b = quantile(out.data1$sims.list$NEC,0.5, names = F)
beta.b = quantile(out.data1$sims.list$beta,0.5, names = F)
top.b = quantile(out.data1$sims.list$top,0.5, names = F)
lines(seq(min.x, NEC.b), rep(quantile(out.data1$sims.list$top,0.5, names = F), length(seq(min.x, NEC.b))))  #this is the NEC line
lines(NEC.b:max.x, top.b*exp(-beta.b*(NEC.b:max.x-NEC.b)))  #this is the decay line
abline(v=NEC.b, lty=2, col = 2)

########################################################
#DRC -NEC
library(drc)
nec.m1 <- drm((suc/tot)~raw.x, data=data1, fct=NEC.4())
summary(nec.m1)   # e = NEC
confint(nec.m1)
pred.m5 <- predict(nec.m1, newdata = df.x, interval="confidence") 

#plotting
plot(data1$raw.x,(data1$suc/data1$tot),log = 'x',main="NEC - DRC") #first plot
lines(vec.x, pred.m5)
abline(v=coef(nec.m1)[4], lty=2, col = 2)


########################################################################################################
#####################################################################################################
#ECx using dose.p and ED + SE

ec = 10 #change ec values here

#10% from the maximum
inh10.1 <- max(df1$prediction)-((ec/100) * max(df1$prediction))
inh10.2 <- max(df2$prediction)-((ec/100) * max(df2$prediction))
inh10.3 <- max(df3$prediction)-((ec/100) * max(df3$prediction))
inh10.4 <- max(df4$Prediction)-((ec/100) * max(df4$Prediction))
inh10.5 <- max(df5$Prediction)-((ec/100) * max(df5$Prediction))
inh10.6 <- max(df6$Prediction)-((ec/100) * max(df6$Prediction))
inh10.7 <- max(df7$prediction)-((ec/100) * max(df7$prediction))
inh10.8 <- max(df8$prediction)-((ec/100) * max(df8$prediction))
inh10.9 <- max(df$prediction)-((ec/100) * max(df$prediction))  #betabinomial

ECfun = function(x,y){
library(VGAM) #Used in EC estimation
eta <- logit(y) 
beta <- coef(x)[1:2] 
ecx <- (eta - beta[1])/beta[2] 
pd <- -cbind(1, ecx)/beta[2] 
ff = as.matrix(vcov(x)[1:2,1:2])
se <- sqrt(((pd %*% ff )* pd) %*% c(1, 1))
upper = (ecx+se*1.96)
lower = (ecx-se*1.96)
sum = data.frame(ecx, lower, upper)
return(sum)
} #using coef. Dervied from dose.p

ECfun2 = function(x,y){
library(VGAM) #Used in EC estimation
  eta <- logit(y) 
  beta <- fixef(x)[1:2] 
  ecx <- (eta - beta[1])/beta[2] 
  pd <- -cbind(1, ecx)/beta[2] 
  ff = as.matrix(vcov(x)[1:2,1:2])
  se <- sqrt(((pd %*% ff )* pd) %*% c(1, 1))
  upper = (ecx+se*1.96)
  lower = (ecx-se*1.96)
  sum = data.frame(ecx, lower, upper)
  return(sum)
} #using fixef

ECfun2 = function(x,y){
  library(VGAM) #Used in EC estimation
  eta <- logit(y) 
  beta <- fixef$cond(x)[1:2] 
  ecx <- (eta - beta[1])/beta[2] 
  pd <- -cbind(1, ecx)/beta[2] 
  ff = as.matrix(vcov(x)[1:2,1:2])
  se <- sqrt(((pd %*% ff )* pd) %*% c(1, 1))
  upper = (ecx+se*1.96)
  lower = (ecx-se*1.96)
  sum = data.frame(ecx, lower, upper)
  return(sum)
} #using fixef

ECfun3 = function(x,y){
  library(VGAM) #Used in EC estimation
  eta <- logit(y) 
  beta <- fixef(x)$cond[1:2] 
  ecx <- (eta - beta[1])/beta[2] 
  pd <- -cbind(1, ecx)/beta[2] 
  ff = as.matrix(vcov(x)$cond[1:2,1:2])
  se <- sqrt(((pd %*% ff )* pd) %*% c(1, 1))
  upper = (ecx+se*1.96)
  lower = (ecx-se*1.96)
  sum = data.frame(ecx, lower, upper)
  return(sum)
} #using fixef

ec1 = ECfun(md1,inh10.1)
ec2 = ECfun(md2,inh10.2 )
ec3 = ECfun2(md3,inh10.3 )
ec4.b = ED(md4, ec, type = c("relative"), interval = "delta")
ec4 = data.frame(ecx = ec4.b[1], lower = ec4.b[3], upper = ec4.b[4])
ec5.b = ED(md5, ec, type = c("relative"), interval = "delta")
ec5 = data.frame(ecx = ec5.b[1], lower = ec5.b[3], upper = ec5.b[4])
ec6.b = ED(md6, ec, type = c("relative"), interval = "delta")
ec6 = data.frame(ecx = ec6.b[1], lower = ec6.b[3], upper = ec6.b[4])
ec7 = ECfun2(md7,inh10.7 )
ec8 = ECfun2(md8,inh10.8 )
ec9 = c(18.8, NA, NA)  #need to change
ec.beta = ECfun3(md3.beta,inh10.9 )
Nec.bayes = quantile(out.data1$sims.list$NEC,c(0.5, .025, .975), names = F)  #Extract median and credible intervals from gamma distribution
NEC.drc = c(coef(nec.m1)[4], confint(nec.m1)[4], confint(nec.m1)[8])

#thresholds 

ec.all =rbind(ec1, ec2, ec3, ec4, ec5, ec6, ec7, ec8, ec9, Nec.bayes, NEC.drc)
#row.names(ec.all) <- c("Binom GLM", 'Quasi GLM', "GLMM", 'LL.4', 'Binom LL.4', 'Weibull','Bayesian GLM','Bayesian GLMM', 'NOEC/ANOVA', 'Nec.bayes', 'NEC.drc')
model <- c("Binom GLM", 'Quasi GLM', "GLMM", 'LL.4 NLR', 'Binom LL.4 NLR', 'Weibull NLR','Bayesian GLM','Bayesian GLMM', 'NOEC/ANOVA', 'Nec.bayes', 'NEC.drc')
ec.all$model = model
ec.all$model <- factor(ec.all$model,levels = c("Binom GLM", "Quasi GLM", "GLMM", "LL.4 NLR", "Binom LL.4 NLR",'Weibull NLR','Bayesian GLM','Bayesian GLMM','NOEC/ANOVA' , 'Nec.bayes', 'NEC.drc')) #Set levels in order
ec.all

#Plotting
library(ggplot2)
p <- ggplot()
p <- p + geom_point(data = ec.all, aes(x = model, y = ecx, size = 6),    col = ("black")) 
p <- p + geom_errorbar(data = ec.all,aes(x = model, ymax = upper, ymin = lower), width=0.2) #add errorbars
p = p + theme(text = element_text(size=6))
p = p + theme_bw()
p = p + coord_flip()
p


# #Plot to check ECvalues
# par(mfrow = c(1, 1))
# plot(data1$raw.x,(data1$suc / data1$tot),log = 'x',main="Checking ECx") #third plot
# lines(df8$vec.x, df8$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
# lines(df8$vec.x, df8$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur upper 95% CI
# lines(df8$vec.x, df8$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur lower 95% CI
# lines(df8$vec.x, rep(inh10.8,lg.x))
# lines(rep(ec8[1], lg.x), seq(0,1,length=lg.x))
# lines(df8$vec.x, rep(inh10.8,lg.x))
# lines(rep(ec8[2], lg.x), seq(0,1,length=lg.x))
# lines(df8$vec.x, rep(inh10.8,lg.x))
# lines(rep(ec8[3], lg.x), seq(0,1,length=lg.x))

##########################################################################################

# #R squared approximations
# library(BaylorEdPsych) #Pseudo R2
# library(MuMIn) #R2 GLMM
# library(ResourceSelection) #Goodness of fit test
# PseudoR2(md1) 
# 1-(md1$deviance / md1$null) # Same as McFadden above
# r.squaredGLMM(md3)

