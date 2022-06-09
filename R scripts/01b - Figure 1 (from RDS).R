
library(ggridges); library(ggplot2); library(magrittr); library(tidyverse); library(vroom)
library(awtools); library(patchwork)

iucn <- readRDS("IUCNDF.rds")

iucn$NVirion[iucn$NVirion==0] <- NA
iucn$NZoon[iucn$NZoon==0] <- NA

iucn %>%
  filter(!(redlistCategory %in% c("Extinct", "Extinct in the Wild"))) %>%
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

################# Tukey tests

## Tukey test

one.way <- aov(iucn$NVirion ~iucn$redlistCategory)
summary(one.way)

tukey.one.way<-TukeyHSD(one.way)
tukey.one.way
plot(tukey.one.way,las=1)
#Tukey Test Zoon#
one.way2 <- aov(iucn$NZoon ~iucn$redlistCategory)
tukey.one.way2<-TukeyHSD(one.way2)
tukey.one.way2

plot(tukey.one.way2,las=1)
