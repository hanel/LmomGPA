library(data.table)
library(ggplot2)
library(lmomRFA)

source('R/auxiliary_functions/imports.R')

trim <- c(0,0)
cluster.number <- 3

{

mx <- data.table(readRDS('data/DV.rds'))

mx <- mx[which(cluster %in% cluster.number),]

mx.trim <- mx[mx[, .I[which(!(dV %in% n.min(dV, trim[1])) & !(dV %in% n.max(dV, trim[2])))], by = SP_ID]$V1]
mx.trim[, p :=  (rank(dV)-.3)/(length(dV) + .4), by = SP_ID]

MX <- data.table(dcast(mx.trim[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

lmom.atsite <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
at.site.para <- t(apply(lmom.atsite, 1, pelgpa))

dta.fit <- sim(MX, dist = 'gpa')
gumbelplot(dta.fit)
qq(dta.fit)

}

s <- sample(dta.fit, length = 500,  type = 'nonpar')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f, method = 'ggplot', rp = F)

s <- sample(dta.fit, length = 500,  type = 'zero')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f, method = 'base', return.period = c(10, 25, 100, 500))

resid.sim <- function(model_object){
  
  dta <- model_object$data
  para <- model_object$REG
  sf <- model_object$scaling_factor
  
  resi <- suppressMessages(melt(dta))
  resi <- data.table(variable = names(dta), sf = sf, t(para))[resi, on = c('variable')]
  resi <- resi[, resi := 1/k*log(1 + k*((value/sf)/alpha)), by = 'variable'] # coles - sigma = alpha, xi = kappa #######

  print(resi[, c('variable', 'resi')])
}

precip <- sim(nim::precip_max[,-1], dist = 'gev')

resid.sim(precip)
resid.sim(dta.fit)

matplot(as.matrix(precip$data), type = 'h')
matplot(as.matrix(dta.fit$data), type = 'h')

precip$scaling_factor
dta.fit$scaling_factor
