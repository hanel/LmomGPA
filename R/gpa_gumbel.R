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

lmom.atsite.t <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
at.site.para <- t(apply(lmom.atsite.t, 1, pelgpa))

# head(lmom.atsite.t)
# head(at.site.para)
# 
# ratiodiagram(lmom.atsite.t[,3:4])
# 
# temp <- as.regdata(data.frame(name = names(MX),
#                               n = apply(MX, 2, function(x) length(which(!is.na(x)))),
#                               mean = lmom.atsite.t[,1],
#                               t = lmom.atsite.t[,2]/lmom.atsite.t[,1],
#                               t_3 = lmom.atsite.t[,3],
#                               t_4 = lmom.atsite.t[,4]))
# 
# para.reg <- regfit(temp, dist = 'gpa')

dta.fit <- sim(MX, dist = 'gpa', trim = c(0, 0))
gumbelplot(dta.fit)
qq(dta.fit)

}

s <- sample(dta.fit, length = 500,  type = 'nonpar')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f, method = 'ggplot', rp = F)

s <- sample(dta.fit, length = 500,  type = 'zero')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f, method = 'base', return.period = c(10, 25, 100, 500))

resid.sim <- function(model_object){}

model_object <- dta.fit

dta <- model_object$data
para <- model_object$REG
sf <- model_object$scaling_factor

resi <- suppressWarnings(melt(dta))
resi <- data.table(variable = names(dta), sf = sf, t(para))[resi, on = c('variable')]
resi <- resi[, resi := 1/k*log(1 + k*((value/sf)/alpha)), by = 'variable']


View(resi)
