
library(ggridges); library(ggplot2); library(magrittr); library(tidyverse); library(vroom)
library(awtools); library(patchwork)

iucn <- read_csv("~/Github/dangerzone/Data/MammalsRedlistVirion_Mar72022.csv")
iucn %<>% filter(!(str_detect(Host, " x ")))
iucn %<>% filter(!(redlistCategory %in% c("Extinct", "Extinct in the Wild")))

nrow(iucn)
table(iucn$redlistCategory)
nrow(iucn[iucn$NVirion > 0,])
table(iucn$redlistCategory[iucn$NVirion > 0])
table(iucn$redlistCategory[iucn$NVirion > 0])/table(iucn$redlistCategory)

table(iucn$redlistCategory[iucn$NZoon>0])
table(iucn$redlistCategory[iucn$NZoon>0])/table(iucn$redlistCategory[iucn$NVirion>0])

#### Make graph, do ANOVA

iucn$NVirion[iucn$NVirion==0] <- NA
iucn$NZoon[iucn$NZoon==0] <- NA

iucn %>%
  mutate(redlistCategory = factor(redlistCategory, levels = rev(c('Data Deficient',
                                                              'Least Concern',
                                                              'Near Threatened',
                                                              'Vulnerable',
                                                              'Endangered',
                                                              'Critically Endangered')))) %>%   filter(!(redlistCategory %in% c("Extinct", "Extinct in the Wild"))) %>%
  pivot_longer(c("NVirion", "NZoon"), names_to = "Level", values_to = "Viruses") %>%
  mutate(Level = recode(Level, !!!c("NVirion"="Total viral richness",
                                    "NZoon" = "Zoonotic viral richness"))) %>%
  ggplot(aes(y = redlistCategory, x = Viruses, fill = redlistCategory)) +
  theme_bw() +
  geom_density_ridges(alpha = 0.5,
                      jittered_points = TRUE,
                      point_alpha=1,
                      point_shape=21) +
  scale_x_continuous(trans = "log10") +
  guides(fill = "none", color = "none") +
  scale_fill_cyclical(values = rev(awtools::a_palette[1:6])) +
  xlab("") + ylab("") + facet_wrap( ~ Level)

################# Simple anova

summary(aov(logVirion ~ redlistCategory, data = iucn[!is.na(iucn$NVirion),]))
summary(aov(logZoon ~ redlistCategory, data = iucn[!is.na(iucn$NVirion),]))

