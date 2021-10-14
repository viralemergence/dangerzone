#path analysis# 

library(lavaan)
library(semPlot)

iucn %<>%
  select(Host, redlistCategory) %<>%
  group_by(Host)

Citations %>%
  select(Host, Citations) %>%
  group_by(Host) %>%
  summarize(Host, Citations) -> citedf
iucn %<>% left_join(citedf)

Path.2
model2 <-'
Citations ~ Endangered + DataDeficient
NVirion ~ Citations + Endangered + Decreasing
NZoon ~ Citations + Endangered + Decreasing
Endangered ~~ DataDeficient
Decreasing ~~ Endangered
NZoon ~~ NVirion'

fit <- cfa(model2, data = FixedPath2)
summary(fit, fit.measures = TRUE, standardized=T,rsquare=T)
semPaths(fit, 'std', layout = 'circle2', curvePivot = TRUE, exoCov = FALSE, residuals = FALSE, 
         nCharNodes = 0, label.cex = 1.3, label.norm = FALSE, label.scale = FALSE, 
         node.width = 2.5, node.height = 2.5, edge.label.cex = 1.1, fade = FALSE)