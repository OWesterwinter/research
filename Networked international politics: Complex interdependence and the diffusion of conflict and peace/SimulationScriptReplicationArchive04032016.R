# Script Monte Carlo simulations JPR article
# "Networked international politics: Complex interdependence
# and the diffusion of conflict and peace" (Dorussen, Gartzke and Westerwinter
##################################################################################

# Set working directory

setwd("C:\\my_working_directory")

# Install and load required libraries

install.packages("sna")
install.packages("ergm")
install.packages("xtable")

library(sna)
library(ergm)
library(xtable)

# Set seed for exact replication of results

set.seed(1243)

# Number of iterations

nsim <- 500

# Vectors in which to store results

full.p.ex <- numeric(nsim)
full.p.end <- numeric(nsim)
omit1.p.ex <- numeric(nsim)
omit1.p.end <- numeric(nsim)
omit2.p.ex <- numeric(nsim)
omit2.p.end <- numeric(nsim)

full.coef.ex <- numeric(nsim)
omit2.coef.ex <- numeric(nsim)

full.coef.se.ex <- numeric(nsim)
omit2.coef.se.ex <- numeric(nsim)

# Simulation

for(i in 1:nsim){

# Create exogenous network

netex <- symmetrize(matrix(rbinom((192*192),1,0.15),nrow=192),rule="strong")

# Compute degree for exogenous network

degex <- degree(netex,gmode="graph",cmode="freeman")

# Create matrix with higher of dyadic degrees

highdegex<-matrix(NA,nrow=192,ncol=192) 

for(j in 1:length(degex)){
highdegex[j,]<-pmax(degex[j],degex)
}

# Rename matrix for model estimation

dyadY1 <- highdegex

# Create exogenous covariate that correlates with higher of
# dyadic degrees

dyadY2<-dyadY1+matrix(rnorm(192*192,sd=1.5),192,192)

# Generate dependent network from distribution with specified parameters

net <- network(192,density=0.01,directed=FALSE)

sim <- simulate(net~edges+edgecov(dyadY1)+kstar(2)+kstar(3)+gwesp(alpha=0,fixed=TRUE), 
coef=c(-5.5,0.45,0.45,-0.25,-0.25),basis=net,
control=control.simulate(MCMC.burnin=5000,MCMC.interval=100))

# Create variable that correlates with dependent network degree

x <- degree(sim,gmode="graph",cmode="freeman")

dyadX <- as.matrix(dist(x,upper=TRUE))+matrix(rnorm(192*192,sd=0.5),192,192)

# Estimate three models

# Full model

full_est <- ergm(sim ~ edges+edgecov(dyadY2)+edgecov(dyadY1)+edgecov(dyadX)+kstar(2)+kstar(3)
+gwesp(alpha=0,fixed=TRUE),estimate=c("MPLE"))

summary(full_est)

# Logit model without exogenous network variable

omit_est1 <- ergm(sim ~ edges+edgecov(dyadY2)+edgecov(dyadX),estimate=c("MPLE"))

summary(omit_est1)

# Logit model with exogenous network variable

omit_est2 <- ergm(sim ~ edges+edgecov(dyadY2)+edgecov(dyadY1)+edgecov(dyadX),estimate=c("MPLE"))

summary(omit_est2)

# Save p-values

full.p.ex[i] <- summary(full_est)$coefs[2,4]

full.p.end[i] <- summary(full_est)$coefs[4,4]

omit1.p.ex[i] <- summary(omit_est1)$coefs[2,4]

omit1.p.end[i] <- summary(omit_est1)$coefs[3,4]

omit2.p.ex[i] <- summary(omit_est2)$coefs[2,4]

omit2.p.end[i] <- summary(omit_est2)$coefs[4,4]

# Save coefficient of exogenous network variable in full model and omit model 2

full.coef.ex[i] <- summary(full_est)$coefs[3,1]

omit2.coef.ex[i] <- summary(omit_est2)$coefs[3,1]

# Save coefficient standard error of exogenous network variable in full and omit model 2

full.coef.se.ex[i] <- summary(full_est)$coefs[3,2]

omit2.coef.se.ex[i] <- summary(omit_est2)$coefs[3,2]

# Count iterations

print(i)

}

# Save results

save(list=c("full.p.ex","full.p.end","omit1.p.ex","omit1.p.end","omit2.p.ex","omit2.p.end",
"full.coef.ex","omit2.coef.ex"),file="SimulationResultsTest03032016.RData")

# Type one error rate (two-tailed p-value) at select p-values

p.vals <- seq(0.01,0.2,by=0.01)

# Calculate type I errors

full.error.ex <- numeric(length(p.vals))
full.error.end <- numeric(length(p.vals))
omit1.error.ex <- numeric(length(p.vals))
omit1.error.end <- numeric(length(p.vals))
omit2.error.ex <- numeric(length(p.vals))
omit2.error.end <- numeric(length(p.vals))

for(i in 1:length(p.vals)){
full.error.ex[i] <- mean(full.p.ex < p.vals[i])
full.error.end[i] <- mean(full.p.end < p.vals[i])
omit1.error.ex[i] <- mean(omit1.p.ex < p.vals[i])
omit1.error.end[i] <- mean(omit1.p.end < p.vals[i])
omit2.error.ex[i] <- mean(omit2.p.ex < p.vals[i])
omit2.error.end[i] <- mean(omit2.p.end < p.vals[i])
}

# Plot results for covariate correlated with exogenous network degree

# Save plot as pdf

pdf("TypeIErrorSimulationExogenousCovariate03032016.pdf",
height=10,width=10)

par(las=1,mar=c(5,5,4,2),cex.lab=1.25,cex.axis=1)
plot(p.vals,full.error.ex,type="l",ylab="P(type I error)",xlab="",ylim=c(0,1))
title("Type I error rate covariate correlated with \n exogenous network variable", line = 1)
abline(h=seq(0,1,by=.1),lty=3,col="grey70")
abline(v=seq(.025,2,by=.025),lty=3,col="grey70")
lines(p.vals,omit1.error.ex,lwd=1.5,col="blue")
lines(p.vals,omit2.error.ex,lwd=1.5,col="green")
points(p.vals,omit1.error.ex,lwd=1.5,col="blue",pch=4,cex=.85)
points(p.vals,full.error.ex,lwd=1.5,pch=1,cex=.85)
points(p.vals,omit2.error.ex,lwd=1.5,col="green",pch=6,cex=.85)
title(xlab="Two-tailed p-value",line=2.25)

legend("right",
legend=c("Logit with exogenous network variable","Logit without exogenous network variable","Full model")[c(3,2,1)],
col=c("green","blue","black")[c(3,2,1)],pch=c(6,4,1)[c(3,2,1)],bg="white")

dev.off()

# Plot results for covariate correlated with endogenous network interdependencies

# Save plot as pdf

pdf("TypeIErrorSimulationEndogenousCovariate03032016.pdf",
height=10,width=10)

par(las=1,mar=c(5,5,4,2),cex.lab=1.25,cex.axis=1)
plot(p.vals,full.error.end,type="l",ylab="P(type I error)",xlab="",ylim=c(0,1))
title("Type I error rate covariate correlated with \n endogenous network interdependence", line = 1)
abline(h=seq(0,1,by=.1),lty=3,col="grey70")
abline(v=seq(.025,2,by=.025),lty=3,col="grey70")
lines(p.vals,omit1.error.end,lwd=1.5,col="blue")
lines(p.vals,omit2.error.end,lwd=1.5,col="green")
points(p.vals,omit1.error.end,lwd=1.5,col="blue",pch=4,cex=.85)
points(p.vals,full.error.end,lwd=1.5,pch=1,cex=.85)
points(p.vals,omit2.error.end,lwd=1.5,col="green",pch=6,cex=.85)
title(xlab="Two-tailed p-value",line=2.25)

legend("topright",
legend=c("Logit with exogenous network variable","Logit without exogenous network variable","Full model")[c(3,2,1)],
col=c("green","blue","black")[c(3,2,1)],pch=c(6,4,1)[c(3,2,1)],bg="white")

dev.off()

# Calculate average estimated beta for exogenous network variable in ergm and logit model

full.coef.ex.avg <- (sum(full.coef.ex))/nsim

omit2.coef.ex.avg <- (sum(omit2.coef.ex))/nsim

# Calculate average bias of estimate for exogenous network variable in ergm and logit model
# (true beta = 0.45)

full.coef.ex.bias <- full.coef.ex.avg-0.45

omit2.coef.ex.bias <- omit2.coef.ex.avg-0.45

# Calculate root mean squared error of ergm and logit model (true beta = 0.45)

full.coef.ex.rmse <- sqrt((sum((full.coef.ex-0.45)^2))/nsim)

omit2.coef.ex.rmse <- sqrt((sum((omit2.coef.ex-0.45)^2))/nsim)

# Create table to summarize bias analysis

full.coef.ex.summary <- c(full.coef.ex.avg,full.coef.ex.bias,full.coef.ex.rmse)

omit2.coef.ex.summary <- c(omit2.coef.ex.avg,omit2.coef.ex.bias,omit2.coef.ex.rmse)

df <- data.frame(full.coef.ex.summary,omit2.coef.ex.summary)

colnames(df) <- c("ERGM","Logit")

rownames(df) <- c("Avg. Beta","Bias","RMSE")

print(df)

# Save table

write.table(df,file="BiasComparisonTable05022016.csv",sep=",",
col.names=c("ERGM","Logit"),row.names=c("Avg. Beta","Bias","RMSE"))

##############################################################
# End of script