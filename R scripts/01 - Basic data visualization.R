library(ggridges); library(ggplot2); library(magrittr); library(tidyverse); library(vroom)

virion <- vroom("~/Github/virion/Virion/virion.csv.gz"); vir <- virion
iucnraw <- read_csv("C:/Users/cjcar/Downloads/IUCNCLOV.csv"); iucn <- iucnraw

vir %<>% filter(HostNCBIResolved = TRUE,
                ICTVRatified = TRUE)

vir %>%
  select(Host, Virus) %>%
  group_by(Host) %>%
  summarize(NVirion = n_distinct(Virus)) -> virdf

iucn %<>% left_join(virdf)

ggplot(iucn, aes(y = redlistCategory, x = NVirion, fill = redlistCategory)) + 
  geom_density_ridges() + 
  scale_x_continuous(trans = "log10")

vir %>% filter(Host == 'homo sapiens') %>%
  select(Virus) %>% distinct() %>% pull(Virus) -> zoonoses

vir %>% filter(Virus %in% zoonoses) %>%
  select(Host, Virus) %>%
  group_by(Host) %>%
  summarize(NZoon = n_distinct(Virus)) -> zoon

iucn %<>% left_join(zoon)

iucn %<>% mutate(redlistCategory = factor(redlistCategory,
                                          levels = rev(c("Data Deficient",
                                                         "Least Concern",
                                                         "Near Threatened",
                                                         "Vulnerable",
                                                         "Endangered",
                                                         "Critically Endangered",
                                                         "Extinct",
                                                         "Extinct in the Wild"))))

ggplot(iucn, aes(y = redlistCategory, x = NZoon, fill = redlistCategory)) + 
  geom_density_ridges() + 
  scale_x_continuous(trans = "log10")

## Tukey test

one.way <- aov(iucn$NVirion ~iucn$redlistCategory)
summary(one.way)

tukey.one.way<-TukeyHSD(one.way)
tukey.one.way
plot(tukey.one.way,las=1)

################# VIRION

#Zoonotic#
virion %>% filter(Host == 'homo sapiens') %>%
  select(Virus) %>% distinct()  %>% pull(Virus) -> zoonoses
virion %>% filter(Virus %in% zoonoses) %>%
  select(Host, Virus) %>%
  group_by(Host) %>%
  summarize(NZoon = n_distinct(Virus)) -> zoon

iucn <- iucnraw 
iucn %<>% left_join(zoon)
iucn %<>% mutate(redlistCategory = factor(redlistCategory,
                                          levels = rev(c("Data Deficient",
                                                         "Least Concern",
                                                         "Near Threatened",
                                                         "Vulnerable",
                                                         "Endangered",
                                                         "Critically Endangered",
                                                         "Extinct",
                                                         "Extinct in the Wild"))))
ggplot(iucn, aes(y = redlistCategory, x = NZoon, fill = redlistCategory)) + 
  geom_density_ridges() + 
  scale_x_continuous(trans = "log10")

#Tukey Test Zoon#
one.way2 <- aov(iucn$NZoon ~iucn$redlistCategory)
tukey.one.way2<-TukeyHSD(one.way2)
tukey.one.way2

plot(tukey.one.way2,las=1)
