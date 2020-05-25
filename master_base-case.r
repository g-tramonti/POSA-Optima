# Optima model, modified for Probabilistic One-way Sensitivity Analysis
#Target Variable: cost of OncotypeDX. See draw_testDX.R for details

rm(list=ls(all=TRUE))
library(MASS)
library(gtools)

#set Global variables
seed <- 10         # set seed for random number generator
S <- 7             # number of states 
disc.b <- 0.035    # discount rate for benefits
disc.c <- 0.035    # discount rate for costs
Nsim <- 10000    # Number of simulations   
lambda <- 20000  #threshold


#load functions
source("draw.R")
source("model.R")
source("draw_testDX.R")
source("draw_testMP.R")
source("draw_testPS.R")
source("draw_testROR_P60.R")

POSA <- function(){

# draw parameters Nsim times
draw()

##-------- MAIN ANALYSIS -----------------

## MCsim for model with chemo for all
# sample from model "Nsim" times
# replace pRec with recurrence probability vector

sim.all <- array(NA,c(Nsim,2))
for (i in 1:Nsim) {
    sim.all[i,]     <- model(i,pRec.all,1,0,)
}
costs.all <- sim.all[,2]
QALYs.all <- sim.all[,1]

## MCsim for model with chemo guided by Oncotype DX
draw.testDX()   # define test-specific parameters
#ctest <- rep(0,Nsim)
sim.DX.high <- array(c(NA,NA),c(Nsim,2))
sim.DX.low  <- array(c(NA,NA),c(Nsim,2))
for (i in 1:Nsim) {
    sim.DX.high[i,] <- model(i,pRec.DXhigh.chemo,1,1,ctestDX)
    sim.DX.low[i,]  <- model(i,pRec.DXlow,0,1, ctestDX)
}
costs.DX <- sim.DX.high[,2]*propHigh + sim.DX.low[,2]*propLow
QALYs.DX <- sim.DX.high[,1]*propHigh + sim.DX.low[,1]*propLow

## MCsim for model with chemo guided by MammaPrint

draw.testMP()   # define test-specific parameters

sim.MP.high <- array(c(NA,NA),c(Nsim,2))
sim.MP.low  <- array(c(NA,NA),c(Nsim,2))
for (i in 1:Nsim) {
  sim.MP.high[i,] <- model(i,pRec.MPhigh.chemo,1,1,ctestMP)
  sim.MP.low[i,]  <- model(i,pRec.MPlow,0,1,ctestMP)
}
costs.MP <- sim.MP.high[,2]*propHighMP + sim.MP.low[,2]*propLowMP
QALYs.MP <- sim.MP.high[,1]*propHighMP + sim.MP.low[,1]*propLowMP


## MCsim for model with chemo guided by Prosignia Subtype

draw.testPS()   # define test-specific parameters

sim.PS.high <- array(c(NA,NA),c(Nsim,2))
sim.PS.low  <- array(c(NA,NA),c(Nsim,2))
for (i in 1:Nsim) {
  sim.PS.high[i,] <- model(i,pRec.PShigh.chemo,1,1,ctestPS)
  sim.PS.low[i,]  <- model(i,pRec.PSlow,0,1,ctestPS)
}
costs.PS <- sim.PS.high[,2]*propHighPS + sim.PS.low[,2]*propLowPS
QALYs.PS <- sim.PS.high[,1]*propHighPS + sim.PS.low[,1]*propLowPS


## MCsim for model with chemo guided by ROR_PT60       #Added

draw.testROR_P60()   # define test-specific parameters

sim.ROR_P60.high <- array(c(NA,NA),c(Nsim,2))
sim.ROR_P60.low  <- array(c(NA,NA),c(Nsim,2))
for (i in 1:Nsim) {
 sim.ROR_P60.high[i,] <- model(i,pRec.ROR_P60high.chemo,1,1,ctestROR_P60)
 sim.ROR_P60.low[i,]  <- model(i,pRec.ROR_P60low,0,1,ctestROR_P60)
}
costs.ROR_P60 <- sim.ROR_P60.high[,2]*propHighROR_P60 + sim.ROR_P60.low[,2]*propLowROR_P60
QALYs.ROR_P60 <- sim.ROR_P60.high[,1]*propHighROR_P60 + sim.ROR_P60.low[,1]*propLowROR_P60



#### Net Benefit analysis  ###

NB.all      = QALYs.all*lambda - costs.all #Net benefits
NB.DX       = QALYs.DX*lambda - costs.DX
NB.MP       = QALYs.MP*lambda - costs.MP
NB.PS       = QALYs.PS*lambda - costs.PS  
NB.ROR_P60  = QALYs.ROR_P60*lambda - costs.ROR_P60   

INMB.DX <- NB.DX - NB.all #Incremental Net Benefits over Chemo for all strategy
INMB.MP <- NB.MP - NB.all
INMB.PS <- NB.PS - NB.all
INMB.ROR_P60 <- NB.ROR_P60 - NB.all

m.NBall<-mean(NB.all) # Expected Net Benefits
m.NBDX<-mean(NB.DX)
m.NBMP<-mean(NB.MP)
m.NBPS<-mean(NB.PS)  
m.NBROR<-mean(NB.ROR_P60)

m.INMB.DX <- mean(INMB.DX) #Expected Incremental Net Benefits
m.INMB.MP <- mean(INMB.MP)
m.INMB.PS <- mean(INMB.PS)
m.INMB.ROR_P60 <- mean(INMB.ROR_P60)

return(c(ctestDX[1],m.NBall,m.NBDX,m.NBMP,m.NBPS,m.NBROR))

}

cNMB <- array(NA,c(11,6))
cent <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)

for (i in 1:11) {
  quant <- cent[i]
  cNMB[i,] <- POSA()
} 
cNMB <- cbind(cent,cNMB)
c.name <- c("Centile","Test cost DX","cNMB Chemo all","cNMB DX","cNMB MP","cNMB PS","cNMB ROR")
colnames(cNMB) <- c.name
write.csv(cNMB,file = "cNMB DX cost.csv")
