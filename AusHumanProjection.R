## Human demographic matrix projections - Australia
## Corey J. A. Bradshaw (corey.bradshaw@flinders.edu.au)
## Flinders University
## Jan 2016

## Remove everything
rm(list = ls())

## libraries
library(boot)

## source
source("matrixOperators.r") 

###############################################
# import data
mort <- read.table("mort.csv", header=T, sep=",")
fert <- read.table("fert.csv", header=T, sep=",")
totpop <- read.table("totpop.csv", header=T, sep=",")
mpop <- read.table("mpop.csv", header=T, sep=",")
fpop <- read.table("fpop.csv", header=T, sep=",")
ltotpop <- dim(totpop)[1]
n14 <- read.table("n14.csv", header=T, sep=",")
Mig.net <- read.table("netMig.csv", header=T, sep=",")
Mig.age <- read.table("Mig.age.avg.exp.csv", header=T, sep=",")
wealth <- read.table("GDP.GPI.csv", header=T, sep=",")
forest <- read.table("forest.csv", header=T, sep=",")

sex.ratio <- mpop[,2:103]/fpop[,2:103]
colnames(sex.ratio) <- colnames(mpop[,2:103])
rownames(sex.ratio) <- mpop[,1] 
sex.ratio2014 <- sex.ratio[44,]
stage.sex.ratio14 <- sex.ratio2014[2:102]

# plots
plot(totpop$yr, totpop$tot, pch=19, type="l")
plot(Mig.net$yr, Mig.net$netMig, type="b", pch=19, cex=0.9, xlab="year", ylab="net immigration")
for.pop <- merge(forest,totpop,by="yr")
for.pop <- for.pop[,1:3]

par(mfrow=c(1,3))
plot(for.pop$tot, for.pop$wvegext, type="b", pch=19, cex=0.8, xlab="population size", ylab="forest extent")
plot(for.pop$yr, for.pop$wvegext, type="b", pch=19, cex=0.8, xlab="year", ylab="forest extent")
plot(for.pop$yr, for.pop$wvegext/for.pop$tot, type="b", pch=19, cex=0.8, xlab="year", ylab="per capita forest extent")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(totpop$tot[-c(ltotpop)], wealth$GDP.US05[-1], type="l", xlab="population size", ylab="GDP (USD 2005)")
plot(totpop$yr[-c(ltotpop)], wealth$GDP.US05[-1]/totpop$tot[-c(ltotpop)], type="l", xlab="year", ylab="per capita GDP (USD 2005)")
plot(totpop$tot[-c(37:ltotpop)], as.vector(na.omit(wealth$GPI.adj.US05[-c(1)])), type="l", xlab="population size", ylab="GPI (USD 2005)")
plot(totpop$yr[-c(37:ltotpop)], as.vector(na.omit(wealth$GPI.adj.US05[-c(1)]))/totpop$tot[-c(37:ltotpop)], type="l", xlab="year", ylab="per capita GPI (USD 2005)")
par(mfrow=c(1,1))

age.vec <- seq(0,100,1)
lage <- length(age.vec)

## construct matrix
stages <- lage
popmat <- matrix(0,nrow=stages,ncol=stages)
colnames(popmat) <- age.vec[1:stages]
rownames(popmat) <- age.vec[1:stages]

## populate matrix
popmat[1, ] <- fert$mf
surv.vec <- 1-mort$Mf
diag(popmat[2:stages, ]) <- surv.vec[-stages]
popmat[stages,stages] <- surv.vec[stages]
popmat.orig <- popmat ## save original matrix

## matrix properties
max.lambda(popmat) ## 1-yr lambda
max.r(popmat) # rate of population change, 1-yr
stable.stage.dist(popmat) ## stable stage distribution
R.val(popmat,stages) # reproductive value
G.val(popmat,stages) # mean generation length

## initial population vector
init.vec <- n14$n14

#################
## project
## set time limit for projection in 1-yr increments
yr.now <- 2014
#************************
yr.end <- 2100 # set projection end date
#************************
t <- (yr.end - yr.now)

## linear fertility trends
int.vec <- 1:t

tot.F <- sum(popmat.orig[1,])
popmat <- popmat.orig


#****************************************************
## Choose scenarios
## Fertility change
F.scen0 <- tot.F ## no fertility change
F.scen1 <- 0.92 # Macao total lowest fertility
F.scen2 <- 1 # worldwide 1-child policy
F.scen3 <- 2.27 # gives lambda = 1
F.scen4 <- 2.03 # UN expectation for 2100
F.scen5 <- 1.54 # doubling of rate of decline compared to UN expectation
F.scen6 <- 2.00 #
#************************
end.fert <- 2100
F.scen.choose <- F.scen0
#************************

#************************
# Net migration
#net.migF <- round(mean(Mig.net$netMigF), 0) # average 2004-2013 (females)
net.migF <- 20000/2

# add migration?
#Mig.scen <- 0 # no net migration
Mig.scen <- 1 # mean annual (fixed) net migration
#Mig.scen <- 2 # random draw of female migrants from data 2004-2013
#Mig.scen <- 3 # constant proportion (0.01) of total population
#Mig.scen <- 4 # double immigration rate as linear increase from 2014
#************************

## weighted average PP emissions of immigrants
imm.ghg.pp.avg <- 5.3468

## Migration change vector (doubles to )
M.mult <- rep(1, t)
M.targ <- 2 # eventual multiplier reached by end.mig
end.mig <- 2100
mt <- (end.mig - yr.now)
M.mult.vec <- pmax(M.mult, (1 - (1 - M.targ)*int.vec/mt)/1)


#****************************************************************************
## Unwanted pregancy aversion
preg.aversion <- 0 # 1 if avoid number of unwanted pregancies; 0 if not avoid
# 2008: 208 million pregnancies; 86 million unintended
# 33 million unplanned births; 41 million abortions; 11 million miscarriages
prop.unwanted <- 33/208 # worldwide proportion of potentially avoidable births
#****************************************************************************
  
## fertilty-change vector
F.scen.ch <- ifelse(F.scen.choose == tot.F, F.scen.choose, F.scen.choose/2)
F.mult <- rep(F.scen.ch/tot.F, t)
ft <- (end.fert - yr.now)
F.mult.vec <- pmax(F.mult, (tot.F - (tot.F - F.scen.ch)*int.vec/ft)/tot.F)

#****************************************************************************
## Stepped fertility reduction (as per reviewer 1's recommended scenario)
## 2.0 child/female by 2025, 1.5 child/female by 2035, and 1 child/female by 2045
stepped <- 0 # 1 if stepped function invoked
fts <- (2045 - yr.now); ints.vec <- 1:fts; Fs.targ <- 1.0*0.5; mults <- Fs.targ/tot.F; F.mults <- rep(mults, fts)
step.fert.mult <- pmax(F.mults, (tot.F - (tot.F - Fs.targ)*ints.vec/fts)/tot.F)
step.fert.mult.vec <- c(step.fert.mult,rep(step.fert.mult[length(step.fert.mult)], ft-fts))
#****************************************************************************
  
#****************************************************************************
## Stepped fertility reduction 2nd scenario
## 2.0 child/female by 2020
stepped2 <- 0 # 1 if stepped function invoked
fts2 <- (2020 - yr.now); ints2.vec <- 1:fts2; Fs.targ2 <- 2.0*0.5; mults2 <- Fs.targ2/tot.F; F.mults2 <- rep(mults2, fts2)
step.fert.mult2 <- pmax(F.mults2, (tot.F - (tot.F - Fs.targ2)*ints2.vec/fts2)/tot.F)
step.fert.mult2.vec <- c(step.fert.mult2,rep(step.fert.mult2[length(step.fert.mult2)], ft-fts2))
#****************************************************************************


## Age at primiparity changes (alpha)
A.scen0 <- 1
A.scen1 <- 0.5 # amount of fertility redistributed from 14:24 to 25:49
#************************
end.alpha <- 2100
A.scen.choose <- A.scen1
#************************

## Non-juvenile (6:oldest) survival change (S)
D.scen0 <- 1 ## no survival change
D.scen1 <- 0.50 # 50 % reduction in stage-specific death rate
D.scen2 <- 0.75 # 25 % reduction in stage-specific death rate
#************************
end.death <- 2100
D.scen.choose <- D.scen0
#************************
## need to create using UN life-expectancy projections (back-calculated to death rates)

## Juvenile survival change (0-5 yrs)
J.scen0 <- 1 ## no survival change
J.scen1 <- 0.50 # 50 % reduction in juvenile death rate
J.scen2 <- 1.50 # 50 % increase in juvenile death rates (e.g., famine from CC)
#************************
end.juvd <- 2100
J.scen.choose <- J.scen0
#************************

## Overall survival changes
S.scen0 <- 1 ## no survival change overall
S.scen1 <- 3.2 # mortality reduces 3.2% across all age classes
S.scen.choose <- S.scen1

## add a catastrophic mortality event
# spread over 5 years
# implemented mid-projection
# equal likelihood of taking any age (/2 for females only)
tot.pop14 <- (totpop[ltotpop,2])
no.toll <- 0
nucl.war.toll <- 1e+9/2 #1,000,000,000
firstsecwars.toll <- 1.31e+8/2 #131,000,000
firstsecwars.prop.toll <- (firstsecwars.toll*2/2500000000*tot.pop14)/2 # deaths proportional to 2.5 b alive at end WWII
superepidemic.toll <- 2e+9/2 #2,000,000,000
killmost.toll <- 6e+9/2 #6,500,000,000
#************************
Cat.scen <- no.toll
#************************

## change to fertility/survival if war/pandemic invoked
war.fert.mult <- 2 # fertility doubles following war/pandemic (subsequently increases linearly until 2013 values thereafter)
war.surv.mult <- 2 # mortality doubles following war (subsequently increases linearly until 2013 values thereafter)

yr.vec <- seq(yr.now,yr.end)

if (Cat.scen != killmost.toll) {  
  eyr.cat <- yr.vec[round(t/2)+5] # year + 1 end of catastrophe
}
if (Cat.scen == killmost.toll) {  
  eyr.cat <- yr.vec[round(t/3)+5] # year + 1 end of catastrophe
}

ws.dt <- (end.death - eyr.cat)
mid.int.vec <- 1:ws.dt
war.surv.mult.vc <- rev(pmin(war.surv.mult, 1 - (1 - war.surv.mult)*mid.int.vec/ws.dt))
# add 1s to first part of vector for no change prior to war/disease catastrophe
add.surv.mult <- rep(1,t-ws.dt)
war.surv.mult.vec <- c(add.surv.mult,war.surv.mult.vc)
war.surv.mult.vec

war.F.mult.vc <- rev(pmin(war.fert.mult, 1 - (1 - war.fert.mult)*mid.int.vec/ws.dt))
add.fert.mult <- rep(1,t-ws.dt)
war.fert.mult.vec <- c(add.fert.mult,war.F.mult.vc)
war.fert.mult.vec

## fertilty-change vector
F.scen.ch <- ifelse(F.scen.choose == tot.F, F.scen.choose, F.scen.choose/2)
F.mult <- rep(F.scen.ch/tot.F, t)
ft <- (end.fert - yr.now)
F.mult.vec <- pmax(F.mult, (tot.F - (tot.F - F.scen.ch)*int.vec/ft)/tot.F)

## changes in female age at primparity (alpha)
## done via redistribution of fertility in 15:24 breeding classes amongst remaining (25:49)
A.mult <- rep(0, t)
at <- (end.alpha - yr.now)
A.mult.vec <- A.scen.choose^(1/t)
A.mult.vec

## changes in survival
## expressed in death rate (1 - S)
D.mult <- rep(1, t)
dt <- (end.death - yr.now)
D.mult.vec <- pmax(D.scen.choose, 1 - (1 - D.scen.choose)*int.vec/dt)
D.mult.vec

## changes in juvenile survival (expressed in mortality changes)
J.mult <- rep(1, t)
jt <- (end.juvd - yr.now)
if (J.scen.choose <= 1) {
  J.mult.vec <- pmax(J.scen.choose, 1 - (1 - J.scen.choose)*int.vec/jt)}
if (J.scen.choose > 1) {
  J.mult.vec <- 1 - (1 - J.scen.choose)*int.vec/jt}
J.mult.vec

## set population storage matrices
n.mat <- matrix(0,nrow=stages,ncol=(t+1))
n.mat[,1] <- init.vec

# immigrant offset emissions vector
mig.e <- rep(0,t-1)

## fertility storage vector
fert.tot <- fert.15.24 <- rep(0,t)

## set up projection loop
for (i in 1:t) {
  fert.tot[i] <- sum(popmat[1,])
  fert.15.24[i] <- sum(popmat[1,16:25])
  n.mat[,i+1] <- popmat %*% n.mat[,i] 
  
  ## add net immigration according to average age breakdown 2004-2013
  if (Mig.scen == 1) {
    n.mat[,i+1] <- n.mat[,i+1] + round((Mig.age$pdF.avg * net.migF), 0)
    migF.e <- sum(round((Mig.age$pdF.avg * net.migF), 0)) * imm.ghg.pp.avg
    migM.e <- sum(round((Mig.age$pdM.avg * net.migF), 0)) * imm.ghg.pp.avg
    mig.e[i] <- migF.e + migM.e
  }
  
  ## add net immigration according to average age breakdown 2004-2013 (random draw)
  if (Mig.scen == 2) {
    net.migF.ran <- sample(Mig.net$netMigF, 1, replace=T)
    n.mat[,i+1] <- n.mat[,i+1] + round((Mig.age$pdF.avg * net.migF.ran), 0)
    migF.e <- sum(round((Mig.age$pdF.avg * net.migF.ran), 0)) * imm.ghg.pp.avg
    migM.e <- sum(round((Mig.age$pdM.avg * net.migF.ran), 0)) * imm.ghg.pp.avg
    mig.e[i] <- migF.e + migM.e
  }
  
  ## add immigrants as constant proportion (0.01) of total population
  if (Mig.scen == 3) {
    n.mat[,i+1] <- n.mat[,i+1] + round((Mig.age$pdF.avg * 0.01*sum(n.mat[,i])), 0)
    migF.e <- sum(round((Mig.age$pdF.avg * 0.01*sum(n.mat[,i])), 0)) * imm.ghg.pp.avg
    male.mig1pc <- 0.01 * sum(n.mat[,i] * as.numeric(stage.sex.ratio14))
    migM.e <- sum(round((Mig.age$pdM.avg * male.mig1pc), 0)) * imm.ghg.pp.avg
    mig.e[i] <- migF.e + migM.e
  }
  
  ## double net immigration according to average age breakdown 2004-2013 (random draw)
    # linear increase from first year of projection
  if (Mig.scen == 4) {
    net.migF.ran <- sample(Mig.net$netMigF, 1, replace=T)
    net.migF.ran.mult <- M.mult.vec[i] * net.migF.ran
    net.migM.ran <- sample(Mig.net$netMigM, 1, replace=T)
    net.migM.ran.mult <- M.mult.vec[i] * net.migM.ran
    n.mat[,i+1] <- n.mat[,i+1] + round((Mig.age$pdF.avg * net.migF.ran.mult), 0)
    migF.e <- sum(round((Mig.age$pdF.avg * net.migF.ran.mult), 0)) * imm.ghg.pp.avg
    migM.e <- sum(round((Mig.age$pdM.avg * net.migM.ran.mult), 0)) * imm.ghg.pp.avg
    mig.e[i] <- migF.e + migM.e
  }  
  
  ## invoke catastrophic mortality over 5-year window
  if (Cat.scen != killmost.toll) {
    if (i == round(t/2,0)) {
      if (Cat.scen == firstsecwars.prop.toll) {
        Cat.scen <- firstsecwars.toll*2/2500000000*sum(n.mat[,i])
      }
      prop.death <- Cat.scen/5/sum(n.mat[,i+1])
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death) ## five-year (Cat.scen/5) catastrophic mortality event
    }
    if (i == round(t/2,0)+1) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
    if (i == round(t/2,0)+2) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
    if (i == round(t/2,0)+3) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
    if (i == round(t/2,0)+4) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
  }

  if (Cat.scen == killmost.toll) {
    if (i == round(t/3,0)) {
      prop.death <- Cat.scen/5/sum(n.mat[,i+1])
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death) ## five-year (Cat.scen/5) catastrophic mortality event
    }
    if (i == round(t/3,0)+1) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
    if (i == round(t/3,0)+2) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
    if (i == round(t/3,0)+3) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
    if (i == round(t/3,0)+4) {
      n.mat[,i+1] <- n.mat[,i+1] - (n.mat[,i+1]*prop.death)
    }
  }
  
  ## unwanted pregnancies averted
  if (preg.aversion == 1) {
    n.mat[1,i+1] <- n.mat[1,i+1] * (1-prop.unwanted) # n * 1-proportion unwanted
  }
  
  ## stepped fertility decline
  if (stepped == 1) {
    F.mult.vec <- step.fert.mult.vec
  }
  
  ## stepped fertility decline (2nd scenario)
  if (stepped2 == 1) {
    F.mult.vec <- step.fert.mult2.vec
  }

  ## change fertility/survival following catastrophic mortality event
  if (Cat.scen > 0) {
    popmat[1,] <- popmat.orig[1,]*war.fert.mult.vec[i]
    diag(popmat[7:stages,6:stages]) <- 1 - ((1 - diag(popmat.orig[7:stages,6:stages]))*war.surv.mult.vec[i])
    popmat[stages,stages] <- 1 - ((1 - popmat.orig[stages,stages])*war.surv.mult.vec[i])
    diag(popmat[2:6,1:5]) <- 1 - ((1 - diag(popmat.orig[2:6,1:5]))*war.surv.mult.vec[i])   
  }
  
  popmat[1,] <- popmat.orig[1,]*F.mult.vec[i]

  # invoke 3.2% reduction in mortality   
  if (S.scen.choose == 3.2) {
    surv.vec.new <- 1-(0.968*mort$Mf)
    diag(popmat[2:stages, ]) <- surv.vec.new[-stages]
    popmat[stages,stages] <- surv.vec.new[stages]
  }
  
  #diag(popmat[7:stages,6:stages]) <- 1 - ((1 - diag(popmat.orig[7:stages,6:stages]))*D.mult.vec[i])
  #popmat[stages,stages] <- 1 - ((1 - popmat.orig[stages,stages])*D.mult.vec[i])
  #diag(popmat[2:6,1:5]) <- 1 - ((1 - diag(popmat.orig[2:6,1:5]))*J.mult.vec[i])
  A.alloc <- sum(popmat[1,16:25])*(1-A.mult.vec)
  F1524.wt <- popmat[1,16:25]/sum(popmat[1,16:25])
  popmat[1,16:25] <- popmat[1,16:25] - F1524.wt*A.alloc
  FR.wt <- popmat[1,26:50]/sum(popmat[1,26:50])
  popmat[1,26:50] <- popmat[1,26:50] + FR.wt*A.alloc 
  ## revisit the mult.vecs - might be able to do faster as for A.alloc
}


## final population size
tot.sex.ratio <- 1/(sum(n14$n14)/(tot.pop14))
m.fin.pop <- sum(n.mat[,(t+1)]*stage.sex.ratio14)
f.fin.pop <- (sum(n.mat[,(t+1)]))
fin.pop <- m.fin.pop + f.fin.pop
#fin.pop <- tot.sex.ratio*(sum(n.mat[,(t+1)]))
fin.pop

## end max lambda
max.lambda(popmat) # 1-yr lambda

# max population size reached over projection interval
#n.max <-  tot.sex.ratio*max(colSums(n.mat))
m.max <- max(rowSums(t(n.mat) %*% diag(as.numeric(stage.sex.ratio14))))
f.max <- max(colSums(n.mat))
n.max <- m.max + f.max
n.max

## x change (initial - final)
m.init.sum <- sum(init.vec*stage.sex.ratio14)
f.init.sum <- sum(init.vec)
times.delta <- round(fin.pop/(m.init.sum+f.init.sum), 2)
times.delta

## year projection vector
yrs <- seq(yr.now,yr.end,1)

## > 65 (or 75, for sensitivity test) proportion (of total population)
m.over65 <- colSums(t(t(n.mat[67:stages,]) %*% diag(as.numeric(stage.sex.ratio14[67:stages]))))
f.over65 <- colSums(n.mat[67:stages,])

ms <- round(as.vector(rowSums(t(n.mat) %*% diag(as.numeric(stage.sex.ratio14)))), 0)
fs <- as.vector(colSums(n.mat))

over65 <- (m.over65+f.over65)/(ms+fs)
#over65 <- colSums(n.mat[67:stages,])/colSums(n.mat)
#over65 <- colSums(n.mat[77:stages,])/colSums(n.mat)

## < 15 proportion (of total pop)
m.under15 <- colSums(t(t(n.mat[1:15,]) %*% diag(as.numeric(stage.sex.ratio14[1:15]))))
f.under15 <- colSums(n.mat[1:15,])

#under15 <- colSums(n.mat[1:15,])/colSums(n.mat)
under15 <- (m.under15+f.under15)/(ms+fs)

## dependency costs
child.cost <- 8500
ret.cost <- 4843

# childcare & aged care costs
under15.cost <- (m.under15+f.under15)*child.cost
over65.cost <- (m.over65+f.over65)*ret.cost
depend.cost <- (under15.cost + over65.cost)

## healthcare costs (2008/2009)
fem.hcc <- c(1870,550,1740,2870,2590,2750,3940,6580,9680,11820)
mal.hcc <- c(2150,760,1070,1250,1580,2410,4120,7250,11670,14310)
age.gr <- c("0-4","5-14","15-24","25-34","35-44","45-54","55-64","65-74","75-84","85+")

# dependent
f0.4 <- colSums(n.mat[1:4,])
f5.14 <- colSums(n.mat[5:14,])
f65.74 <- colSums(n.mat[65:74,])
f75.84 <- colSums(n.mat[75:84,])
f85p <- colSums(n.mat[85:stages,])

m0.4 <- colSums(t(t(n.mat[1:4,]) %*% diag(as.numeric(stage.sex.ratio14[1:4]))))
m5.14 <- colSums(t(t(n.mat[5:14,]) %*% diag(as.numeric(stage.sex.ratio14[5:14]))))
m65.74 <- colSums(t(t(n.mat[65:74,]) %*% diag(as.numeric(stage.sex.ratio14[65:74]))))
m75.84 <- colSums(t(t(n.mat[75:84,]) %*% diag(as.numeric(stage.sex.ratio14[75:84]))))
m85p <- colSums(t(t(n.mat[85:stages,]) %*% diag(as.numeric(stage.sex.ratio14[85:stages]))))

depend.fhcost <- ((f0.4)*fem.hcc[1])+((f5.14)*fem.hcc[2])+((f65.74)*fem.hcc[8])+((f75.84)*fem.hcc[9])+((f85p)*fem.hcc[10])
depend.mhcost <- ((m0.4)*mal.hcc[1])+((m5.14)*mal.hcc[2])+((m65.74)*mal.hcc[8])+((m75.84)*mal.hcc[9])+((m85p)*mal.hcc[10])
depend.hcost <- (depend.fhcost + depend.mhcost)

# working health costs (not used)
f15.24 <- colSums(n.mat[15:24,])
f25.34 <- colSums(n.mat[25:34,])
f35.44 <- colSums(n.mat[35:44,])
f45.54 <- colSums(n.mat[45:54,])
f55.64 <- colSums(n.mat[55:64,])


## constant working population gdp
work.gdp <- 66409
m.working <- colSums(t(t(n.mat[16:66,]) %*% diag(as.numeric(stage.sex.ratio14[16:66]))))
f.working <- colSums(n.mat[16:66,])
working.gdp <- (m.working+f.working) * work.gdp


## dependency ratio
# number of people < 15 and > 65 relative to rest
# more older people offset by fewer dependent young people (paper by Ehrlichs)
#dep.ratio <- (colSums(n.mat[1:15,]) + colSums(n.mat[67:stages,])) / colSums(n.mat[16:66,])
dep.ratio <- ((m.under15+f.under15) + (m.over65+f.over65)) / (m.working+f.working)
depend.cost.ratio <- (depend.cost+depend.hcost)/working.gdp


## plots
par(mfrow=c(2,2),yaxt="s")
plot(yrs,as.vector(colSums(n.mat)*tot.sex.ratio),type="l",xlab="year",ylab="N",ylim=c(0,1.02*n.max),xlim=c(yr.now,yr.end))
abline(h=sum(n.mat[,1]*tot.sex.ratio),lty=2)
maxN.sub <- which(colSums(n.mat)==max(colSums(n.mat)))
abline(v=yrs[maxN.sub],lty=2)
title(main=paste("final N = ",round(fin.pop/1e6, 3), " m","; max N = ",round(n.max/1e6, 3), " m",sep=""),
      sub = paste("delta = ", times.delta, "x", sep=""))

plot(yrs,as.vector(over65),type="l",xlab="year",ylab="proportion",ylim=c(0,1))
lines(yrs,as.vector(under15),lty=2)
title(main=paste("prop <15 (dashed) & >65 (solid) yrs", sep=""))

plot(yrs,as.vector(dep.ratio),type="l",xlab="year",ylab="dependency ratio",ylim=c(0,1))
title(main=paste("dependency ratio: initial = ", round(dep.ratio[1],4), "; final = ", round(dep.ratio[t+1],4), sep=""))

plot(yrs,as.vector(depend.cost.ratio),type="l",xlab="year",ylab="dependency cost ratio", ylim = c(0.0,1))
title(main=paste("dependency cost ratio to working GDP: initial = ", round(depend.cost.ratio[1],4), "; final = ", round(depend.cost.ratio[t+1],4), sep=""))

par(mfrow=c(1,1))

## average annual growth 1971-2014
lam.totpop <- totpop$tot[2:ltotpop]/totpop$tot[1:(ltotpop-1)]
mlam7114 <- mean(lam.totpop)
mlam7114
mr7114 <- mean(log(lam.totpop))
mr7114

## if 2014 population grew at same mean rate 1971-2014, there would be x by 2100
npred.2100 <- totpop$tot[ltotpop] * mlam7114^t
npred.2100

## average annual growth 2015-2100
lam.nmat <- as.vector(fs+ms)[2:(t+1)]/as.vector(fs+ms)[1:t]
mean(lam.nmat)

## make new plot with previous data + projection
years <- c(totpop$yr[-ltotpop], yrs)
ns <- c(totpop$tot[-ltotpop], as.vector(ms+fs))
  
plot(years, ns, type="l", xlab="year", ylab="total population")
points(totpop$yr, totpop$tot, pch=19, cex=0.6)
title(main=paste("final N = ",round(fin.pop/1e6, 3), " m","; max N = ",round(n.max/1e6, 3), " m",sep=""), sub = paste("delta = ", times.delta, "x; ", "mn.lam(2014-2100) = ", round(mean(lam.nmat),5), sep=""))
abline(v=2014, lty=2)
abline(h=totpop$tot[ltotpop], lty=2)


## export for plotting
# dependency
 depratio <- as.vector(dep.ratio)
 O65 <- as.vector(over65); U15 <- as.vector(under15)
 outd <- data.frame(yrs,depratio,O65,U15)
 write.table(outd,"depend.avg2immig.csv",sep=",", row.names = F, col.names = TRUE)

# dependency cost ratio
depcostratio <- as.vector(depend.cost.ratio)
outdcost <- data.frame(yrs,depcostratio)
write.table(outdcost,"depend.cost.20Kimmig.310316.csv",sep=",", row.names = F, col.names = TRUE)

# N
 out <- data.frame(years,ns)
 write.table(out,"immig2avg.out.csv",sep=",", row.names = F, col.names = TRUE)




############################
## post-hoc carbon estimates
############################
ghg.pp.avg <- 26.6221 # tonnes CO2-e per person (average 1990-2012, NGGI+LULUCF)
#ghg.pp.avg <- 24.4000 # tonnes CO2-e per person (average 1990-2012, NGGI)

ghg.pp.2050.target <- 5.818883e+008*0.2 # 80% below 2000 level (NGGI+LULUCF)
ghg.pp.2020.target <- 5.818883e+008*0.95 # 5% below 2000 level (NGGI+LULUCF)
ghg.pp.2030.target <- 6.0695726153e+008*(1-0.27) # 27% below 2005 level (NGGI+LULUCF)
 
1 - (ghg.pp.2030.target/(ghg.pp.2020.target/0.95)) # 2030 target expressed in 2000 (not 2005) baseline terms

# 2030
ghg.pp.avg.series <- ns[43:60]*ghg.pp.avg # 2030
out.ghg.pp.avg.series <- data.frame(years[43:60],ghg.pp.avg.series) # 2030
mig.co2e.offset <- data.frame(yrs[2:17],mig.e[1:16]) ## immigrant emissions offset
mig.co2e.offset.sum <- sum(out.ghg.pp.avg.series[3:18,2] - mig.co2e.offset[,2]) ## net CO2 
ghg.pp.avg.series[16] # 2030
ghg.pp.2030 <- ghg.pp.2030.target/ns[60]  # per capita emissions in 2020
ghg.pp.2030

# output values for graphs
mig.co2e.offset.sum/10^9 # Gt cumulative net emissions after removing immigrant emissions prior to immigration
sum(ghg.pp.avg.series[3:18])/10^9 # Gt cumulative emissions

# 2050
ghg.pp.avg.series <- ns[43:80]*ghg.pp.avg # 2050
out.ghg.pp.avg.series <- data.frame(years[43:80],ghg.pp.avg.series) # 2050
mig.co2e.offset <- data.frame(yrs[2:37],mig.e[1:36]) ## immigrant emissions offset
mig.co2e.offset.sum <- sum(out.ghg.pp.avg.series[3:38,2] - mig.co2e.offset[,2]) ## net CO2 
ghg.pp.avg.series[36] # 2050
ghg.pp.2050 <- ghg.pp.2050.target/ns[80] # per capita emissions in 2050
ghg.pp.2050

# output values for graphs
mig.co2e.offset.sum/10^9 # Gt cumulative net emissions after removing immigrant emissions prior to immigration
sum(ghg.pp.avg.series[3:38])/10^9 # Gt cumulative emissions


ghg.pp.2012 <- 26.6221 # total
#ghg.pp.2012 <- 24.4000 # NGGI

# 2030
ghg.pp.2012/ghg.pp.2030
ghg.pp.pyr <- (ghg.pp.2012 - ghg.pp.2030)/(2030-2012)
pp.red.seq <- seq(ghg.pp.2012,ghg.pp.2030,-ghg.pp.pyr)
out.pp.red.seq <- data.frame(years[42:60],pp.red.seq)
out.pp.red.seq

# 2050
ghg.pp.2012/ghg.pp.2050
ghg.pp.pyr <- (ghg.pp.2012 - ghg.pp.2050)/(2050-2012)
pp.red.seq <- seq(ghg.pp.2012,ghg.pp.2050,-ghg.pp.pyr)
out.pp.red.seq <- data.frame(years[42:80],pp.red.seq)
out.pp.red.seq

write.table(out.pp.red.seq,"ghg.pp.red.seq.avgmig.out.csv",sep=",", row.names = F, col.names = TRUE)
