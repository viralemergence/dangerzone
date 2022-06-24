#path analysis#

library(tidyverse)
library(magrittr)
library(lavaan)
library(semPlot)

iucn <- read_csv("~/Github/dangerzone/Data/FixedPathMar72022.csv")
iucn %<>% filter(!(str_detect(Host, " x ")))

iucn %<>% filter(complete.cases(NVirion, NZoon, Endangered, DataDeficient, Decreasing, Citations))
# All values are already log transformed

iucn %<>% rename(Viruses = NVirion,
                 Zoonotic = NZoon,
                 DataDef = DataDeficient)

iucn %<>% filter(Viruses > 0)

model2 <-'
Citations ~ Endangered
DataDef ~ Citations
Viruses ~ Citations + Endangered + Decreasing + DataDef
Zoonotic ~ Citations + Endangered + Decreasing + DataDef
Endangered ~~ DataDef
Decreasing ~~ Endangered
Zoonotic ~~ Viruses'

fit <- cfa(model2, data = iucn)
summary(fit, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit, what = 'est', layout = 'circle2', curvePivot = TRUE, exoCov = FALSE, residuals = FALSE,
         nCharNodes = 0, label.cex = 1, label.norm = FALSE, label.scale = FALSE,
         node.width = 1.7, node.height = 1.7, edge.label.cex = 1.1, fade = FALSE)
