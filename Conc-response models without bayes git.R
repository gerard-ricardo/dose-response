#all models  
hav4 = read.table("https://pastebin.com/raw/v7Uy4iXs", header= TRUE,dec=",")
hav4

#options(scipen = 999)  # turn off scientific notation

# Load packages (if not installed Tools->Install Packages)
library(lattice)
library(MASS) 
library(lme4) #GLM
library(drc) #NLR
library(RVAideMemoire) #GLMM overdispersion test
library(VGAM) #Used in EC estimation
library(BaylorEdPsych) #Pseudo R2
library(MuMIn) #R2 GLMM
library(ResourceSelection) #Goodness of fit test
library(ggplot2) #GGplot
library(reshape2)

##########################################
#Data organising
hav4$raw.x <- as.numeric(as.character(hav4$raw.x))
hav4$suc <- as.numeric(as.character(hav4$suc))
hav4$tot <- as.numeric(as.character(hav4$tot))

#setting predicition
min.x <- min(hav4$raw.x)
max.x <- max(hav4$raw.x)
lg.x <- 1000
df.x <- data.frame(raw.x = seq(min.x, max.x, length = lg.x)) #setting up  new  data frame (df) defining log.x values to run 
vec.x =df.x[,1]

#setting space plots
par(mfrow = c(3, 2)) 


########################################################################################
#1) Binomial GLM 
md1 <- glm(cbind(suc,(tot - suc)) ~ raw.x, family = binomial (link = logit),data = hav4)
summary(md1)
sum(residuals(md1, type = "deviance")^2)/md1$df.residual # >1.5, overdispersed


GLMfun = function(i){
  mm <- model.matrix(~raw.x, df.x)  # build model matrix 
  eta <- mm %*% coef(i)
  prediction  <- as.vector(exp(eta) / (1 + exp(eta)))
  se    <- sqrt(diag(mm %*% vcov(i) %*% t(mm)))
  upper  <- exp(eta + 1.96 *se) /(1 + exp(eta  + 1.96 *se))
  lower  <- exp(eta - 1.96 *se) /(1 + exp(eta  - 1.96 *se))
  df = data.frame(vec.x, prediction, lower,upper)
  return(df)
}
df1=GLMfun(md1)

plot(hav4$raw.x,(hav4$suc / hav4$tot),log = 'x', main="Binomial GLM") #first plot
lines(df1$vec.x, df1$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df1$vec.x, df1$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df1$vec.x, df1$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

########################################################################################
#2) Run a quasibinomial GLM
md2 <- glm(cbind(suc,(tot - suc)) ~ raw.x, family = quasibinomial (link = logit),data = hav4)

df2=GLMfun(md2)

plot(hav4$raw.x,(hav4$suc / hav4$tot),log = 'x', main="Quasibinomial GLM") #second plot
lines(df2$vec.x, df2$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df2$vec.x, df2$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df2$vec.x, df2$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

#######################################################################################################
#3) Run a GLMM with a random slope effect (i.e tank / container etc)
hav4$obs <- factor(formatC(1:nrow(hav4), flag="0", width = 3))# first need to make an observation row to soon randomise
md3 <- glmer(cbind(suc,(tot - suc)) ~ raw.x + (1|obs) ,family = binomial (link = logit),data = hav4) 
summary(md3)
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
}
df3 = GLMfun2(md3)

#plot
plot(hav4$raw.x,(hav4$suc / hav4$tot),log = 'x',main="GLMM") #third plot
lines(df3$vec.x, df3$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df3$vec.x, df3$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur upper 95% CI
lines(df3$vec.x, df3$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur lower 95% CI

############################################################################################
# 4) Non-linear regression
md4 <- drm((suc / tot)~ raw.x, fct = LL.4(fixed = c(NA, 0, NA, NA),names = c("Slope", "Lower", "Upper", "ED50")),data = hav4) #constrained to zero

pred.m4 <- predict(md4, newdata = df.x, interval="confidence") #predict 95% CI
df4 = data.frame(vec.x, pred.m4)

plot(hav4$raw.x,(hav4$suc/ hav4$tot),log="x", main="LL.4 Lower=0") #fouth plot
lines(df4$vec.x, df4$Prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df4$vec.x, df4$Lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df4$vec.x, df4$Upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

############################################################################################################
# 5) Non-linear regression (binomial)
md5 <- drm((suc / tot)~ raw.x, fct = LL.4(fixed = c(NA, 0, NA, NA),names = c("Slope", "Lower", "Upper", "ED50")),  weights = tot, type = 'binomial',data = hav4) 

pred.m5 <- predict(md5, newdata = df.x, interval="confidence") #predict 95% CI
df5 = data.frame(vec.x, pred.m5)

plot(hav4$raw.x,(hav4$suc/ hav4$tot),log="x", main="LL.4 Binomial Lower=0") #fouth plot
lines(df5$vec.x, df5$Prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df5$vec.x, df5$Lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df5$vec.x, df5$Upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI

##################################################################################
# 6) Non-linear regression
md6 <- drm((suc / tot)~ raw.x, fct = W1.4(fixed = c(NA, 0, NA, NA),names = c("Slope", "Lower", "Upper", "ED50")),weights = tot, type = 'binomial',data = hav4) 

pred.m6 <- predict(md6, newdata = df.x, interval="confidence") #predict 95% CI
df6 = data.frame(vec.x, pred.m6)

plot(hav4$raw.x,(hav4$suc/ hav4$tot),log="x", main="W1.4 Weibull Lower=0") #fouth plot
lines(df6$vec.x, df6$Prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df6$vec.x, df6$Lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot upper 95% CI
lines(df6$vec.x, df6$Upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #plot lower 95% CI


#EC10 using dose.p and ED + SE

#10% from the maximum
inh10.1 <- max(df1$prediction)-(0.1 * max(df1$prediction))
inh10.2 <- max(df2$prediction)-(0.1 * max(df2$prediction))
inh10.3 <- max(df3$prediction)-(0.1 * max(df3$prediction))
inh10.4 <- max(df4$Prediction)-(0.1 * max(df4$Prediction))
inh10.5 <- max(df5$Prediction)-(0.1 * max(df5$Prediction))
inh10.6 <- max(df6$Prediction)-(0.1 * max(df6$Prediction))

ECfun = function(x,y){
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
} #using coef
ECfun2 = function(x,y){
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

ec1 = ECfun(md1,inh10.1)
ec2 = ECfun(md2,inh10.2 )
ec3 = ECfun2(md3,inh10.3 )
ec4.b = ED(md4, 10, type = c("relative"), interval = "delta")
ec4 = data.frame(ecx = ec4.b[1], lower = ec4.b[3], upper = ec4.b[4])
ec5.b = ED(md5, 10, type = c("relative"), od = T, interval = "delta")
ec5 = data.frame(ecx = ec5.b[1], lower = ec5.b[3], upper = ec5.b[4])
ec6.b = ED(md6, 10, type = c("relative"), interval = "delta")
ec6 = data.frame(ecx = ec6.b[1], lower = ec6.b[3], upper = ec6.b[4])
ec.all =rbind(ec1, ec2, ec3, ec4, ec5, ec6)
row.names(ec.all) <- c("Binom GLM", 'Quasi GLM', "GLMM", 'LL.4', 'Binom LL.4', 'Weibull')
model <- c("Binom GLM", 'Quasi GLM', "GLMM", 'LL.4', 'Binom LL.4', 'Weibull')
ec.all$model = model
ec.all$model <- factor(ec.all$model,
                         levels = c("Binom GLM", "Quasi GLM", "GLMM", "LL.4", "Binom LL.4",'Weibull')) #Set levels in order
ec.all
#Plotting
p <- ggplot()
p <- p + geom_point(data = ec.all, aes(x = model, y = ecx, size = 6),    col = ("black")) 
p <- p + geom_errorbar(data = ec.all,aes(x = model, ymax = upper, ymin = lower), width=0.2) #add errorbars
p = p + theme(text = element_text(size=6))
p = p + theme_bw()
p

#Plot to check ECvalues (using GLMM as an exmaple)
par(mfrow = c(1, 1))
plot(hav4$raw.x,(hav4$suc / hav4$tot),log = 'x',main="Checking ECx") #third plot
lines(df3$vec.x, df3$prediction, type = "l", lwd = 2, col = 2, xaxt = "n", las = 1) #plot model mean line
lines(df3$vec.x, df3$upper, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur upper 95% CI
lines(df3$vec.x, df3$lower, type = "l", lwd = 2,  xaxt = "n", las = 1) #Zuur lower 95% CI
lines(df3$vec.x, rep(inh10.3,lg.x))
lines(rep(ec3[1], lg.x), seq(0,1,length=lg.x))
lines(df3$vec.x, rep(inh10.3,lg.x))
lines(rep(ec3[2], lg.x), seq(0,1,length=lg.x))
lines(df3$vec.x, rep(inh10.3,lg.x))
lines(rep(ec3[3], lg.x), seq(0,1,length=lg.x))

##########################################################################################

#R squared approximations
PseudoR2(md1) #GLMs
1-(md1$deviance / md1$null) # Same as McFadden above
r.squaredGLMM(md3)



#######################################################


