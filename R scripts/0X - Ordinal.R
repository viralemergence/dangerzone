#path analysis#

library(tidyverse)
library(magrittr)
library(ordinal)

iucn <- readRDS("IUCNDF.rds")

iucn %<>% filter(complete.cases(NVirion, NZoon, Endangered, Data_Deficient, Decreasing, Citations))
iucn %<>% mutate(Citations = log(Citations + 1))
iucn %<>% mutate(NVirion = log(NVirion))
iucn %<>% mutate(NZoon = log(NZoon))


iucn %<>% rename(Viruses = NVirion,
                 Zoonotic = NZoon,
                 DataDef = Data_Deficient)
# iucn %<>%
#   select(Host, redlistCategory) %<>%
#   group_by(Host)
#
# Citations %>%
#   select(Host, Citations) %>%
#   group_by(Host) %>%
#   summarize(Host, Citations) -> citedf
# iucn %<>% left_join(citedf)

## A simple cumulative link model:
# fm1 <- clm(rating ~ contact + temp, data=wine)
# summary(fm1)

## A simple cumulative link mixed model:
# fmm1 <- clmm(rating ~ contact + temp + (1|judge), data=wine)
# summary(fmm1)

fm1 <- clm()
