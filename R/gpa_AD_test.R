library(data.table)
library(ggplot2)
library(lmom)

source('R/aux_functions.R')

trim <- c(0,0)
smp <- 1000

mx <- data.table(readRDS('data/DV.rds'))

MX <- as.data.table(dcast(mx[, list(year, SP_ID, dV)], year ~ SP_ID, value.var = 'dV'))
MX <- MX[, year := NULL]

# lmom.atsite.t <- as.data.frame(t(apply(MX, 2, function(x) samlmu(x, trim = trim))))
# gg.MRD(lmom.atsite.t[,3:4]) + ggtitle('via lmom package')
# at.site.pars <- data.table(SP_ID = names(MX), t(apply(lmom.atsite.t, 1, function(x) {pelgpa(x)})))

lmom.atsite.manual <- as.data.frame(t(apply(MX, 2, function(x) moje.lm(x)))) # at site l-momenty (prepsanou fci)
gg.MRD(lmom.atsite.manual[,3:4]) + ggtitle('via moje l-momenty')
at.site.pars.manual <- data.table(SP_ID = names(MX), t(apply(lmom.atsite.manual, 1, function(x) {gpa.para(x)}))) # atsite gpa parametry (prepsanou fci)

# identical(summary(at.site.pars),summary(at.site.pars.manual))

para.mx <- merge(mx[, .(SP_ID, dV)], at.site.pars.manual)
ad.mx <- para.mx[, .(base_ad = AD.test.multidist(val = dV, # vypocet AD statistiky pro deficitni objemy
                                                 location = unique(xi), 
                                                 scale = unique(alpha), 
                                                 shape = unique(k), 
                                                 distr = 'gpa')), 
                 by = SP_ID]

AD.SMP <- lapply(1:smp, function(i) { 
  ad.dV <- para.mx[, .(dV = quagpa(runif(length(dV)), c(unique(xi), unique(alpha), unique(k)))), by = SP_ID] # samplovani GPA hodnot
  ad.para <- dcast(ad.dV[, .(val = gpa.para(moje.lm(dV)), para = c('xi', 'alpha', 'k')), by = SP_ID], SP_ID ~ para, value.var = 'val') # paramety samplu
  ad.res <- merge(ad.dV, ad.para)
  ad.res[, .(smp_ad = AD.test.multidist(val = dV, location = unique(xi), scale = unique(alpha), shape = unique(k), distr = 'gpa')), by = SP_ID] # AD stat. pro samply
})

ad.smp <- do.call(rbind, AD.SMP)

ad <- merge(ad.mx, ad.smp)
res <- unique(ad[, .(p = length(which(smp_ad >= base_ad))/ .N, # vypocty p-hodnoy a vyhodnoceni testu s hl. vyzn. 5%
                     eval = base_ad < quantile(smp_ad, .95)), by = SP_ID]) #######################

length(which(res$eval))

summary(res$p)
summary(res$eval)

gg.MRD(lmom.atsite.manual[which(!res$eval),3:4]) + ggtitle('ty co neprosly testem')

# saveRDS(res, 'Active Docs/filip_nim/data/AD.rds')
#############################################################

library(rgdal)
library(maptools)

pov <- readOGR('data/mezipovodi_133_cluster.shp', 'mezipovodi_133_cluster', verbose = F)

centroid <- data.frame(SP_ID = pov@data$DBCN, coordinates(pov)) 

pov.f <- fortify(pov, region = 'DBCN')
names(pov.f)[which(names(pov.f) == 'id')] <- 'SP_ID'

polygony <- Reduce(function(...) merge(..., all = TRUE, by = 'SP_ID'), list(pov.f, centroid, res))

ad.map <- ggplot(polygony) +
  geom_polygon(aes(x = long, y = lat, group = SP_ID, fill = factor(eval)), color = 'grey25') +
  geom_text(aes(x = X1, y = X2, group = SP_ID,label = round(p, 2)), size = 2.5) +
  scale_fill_manual(values = c('red4','lightsteelblue'), name = '') +
  coord_fixed(ratio = 1) + 
  theme_classic()

ad.map
# htmlwidgets::saveWidget(plotly::as_widget(plotly::ggplotly(ad.map)), file = paste0(getwd(),'/Active Docs/filip_nim/data/ad_map.html'))
