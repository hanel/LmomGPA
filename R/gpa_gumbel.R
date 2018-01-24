library(data.table)
library(ggplot2)
library(lmomRFA)

source('R/auxiliary_functions/imports.R')

trim <- c(1,2)
cluster.number <- 1

{

mx <- data.table(readRDS('data/DV.rds'))

mx <- mx[which(cluster %in% cluster.number),]

mx.trim <- mx[mx[, .I[which(!(dV %in% n.min(dV, trim[1])) & !(dV %in% n.max(dV, trim[2])))], by = SP_ID]$V1]
mx.trim[, p :=  (rank(dV)-.3)/(length(dV) + .4), by = SP_ID]

MX <- data.table(dcast(mx.trim[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

lmom.atsite <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
at.site.para <- t(apply(lmom.atsite, 1, pelgpa))

}

dta.fit <- sim(MX[,2], dist = 'gpa')
gumbelplot(dta.fit)
qq(dta.fit)

geom_qq()

s <- sample(dta.fit, length = 500,  type = 'nonpar')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f, return.period = c(10, 25, 100, 500))

s <- sample(dta.fit, length = 500,  type = 'zero')
f <- fit(s, dta.fit)

growthcurve(dta.fit, f, method = 'ggplot', rp = F)


res <- resid.sim(dta.fit)
res
xxx <- backtodata(res, dta.fit)
xxx
