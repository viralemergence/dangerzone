library(ggridges); library(ggplot2); library(magrittr); library(tidyverse); library(vroom)
library(awtools); library(patchwork)

# virion <- vroom("~/Github/virion/Virion/virion.csv.gz"); vir <- virion
# iucnraw <- read_csv("C:/Users/cjcar/Downloads/IUCNCLOV.csv"); iucn <- iucnraw

vir <-
  virion <-
  "Data/Virion.csv.gz" %>% vroom::vroom()

iucn <-
  iucnraw <-
  read_csv("Data/iucnclov.csv")

vir %<>% filter(HostNCBIResolved = TRUE,
                ICTVRatified = TRUE)

vir %>%
  select(Host, Virus) %>%
  group_by(Host) %>%
  summarize(NVirion = n_distinct(Virus)) -> virdf

iucn %<>% left_join(virdf)

# ggplot(iucn, aes(y = redlistCategory, x = NVirion, fill = redlistCategory)) +
#   geom_density_ridges() +
#   scale_x_continuous(trans = "log10")

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
# iucn %>%
#   filter(!(redlistCategory %in% c("Extinct", "Extinct in the Wild"))) %>%
#   ggplot(aes(y = redlistCategory, x = NVirion, fill = redlistCategory)) +
#   theme_ridges() +
#   geom_density_ridges(alpha = 0.5,
#                       jittered_points = TRUE,
#                       point_alpha=1,
#                       point_shape=21) +
#   scale_x_continuous(trans = "log10") +
#   guides(fill = "none", color = "none") +
#   scale_fill_cyclical(values = rev(awtools::a_palette[1:6])) +
#   xlab("") + ylab("") -> g1

# Zoonotic

# iucn %>%
#   filter(!(redlistCategory %in% c("Extinct", "Extinct in the Wild"))) %>%
#   pivot_longer(c("NVirion", "NZoon"), names_to = "Level", values_to = "Viruses") %>%
#   mutate(Level = recode(Level, !!!c("NVirion"="Total viral richness",
#                                     "NZoon" = "Zoonotic viral richness"))) %>%
#   ggplot(aes(y = redlistCategory, x = Viruses, fill = redlistCategory)) +
#   theme_bw() +
#   geom_density_ridges(alpha = 0.5,
#                       jittered_points = TRUE,
#                       point_alpha=1,
#                       point_shape=21) +
#   scale_x_continuous(trans = "log10") +
#   guides(fill = "none", color = "none") +
#   scale_fill_cyclical(values = rev(awtools::a_palette[1:6])) +
#   xlab("") + ylab("") + facet_wrap( ~ Level)

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
