
# 03_Greg Analysis ####

rm(list = ls())

{

  library(tidyverse); library(data.table); library(magrittr); library(ggregplot); library(cowplot)
  library(geiger);library(ape);library(picante)
  library(tidyverse); library(INLA); library(ggregplot); library(glue); library(fs)

  theme_set(theme_cowplot())

}

iucn <- readRDS("IUCNDF.rds")

iucn %<>% filter(!populationTrend == "")

# source("0_Urban Import.R")

dir_create("Output")

# UrbanDF %<>% filter(!Pseudoabsence)

Resps <- c(

  "NVirion", "NZoon"

)

FamilyList <- rep("nbinomial", 2)

names(FamilyList) <- Resps

FullCovar <- c("hOrder",
               "LogCites",
               "Endangered", "Data_Deficient",
               "populationTrend",
               # "DomesticBinary", "HumanDistance",
               "LogArea", "DietDiversity",
               "LogMass",
               "PC1", "PC2")

PhyloList <-
  IMList <-
  TestDFList <-
  list()

# iucn %<>%
#   filter(!is.na(NVirion)) %>% # Removing all pseudoabsences
#   filter(!is.na(NVirion) & !(LogCites == 0)) %>% # Removing pseudoabsences and citations
#   mutate_at("NZoon", ~ifelse(is.na(.x), 0, .x)) %>%
#   mutate_at("NVirion", ~.x - 1)

iucn %<>%
  mutate_at(c("NZoon", "NVirion"), ~ifelse(is.na(.x), 0, .x))

iucn %<>% mutate(LogRichness = log(NVirion + 1))

q <- 1

for(q in q:length(Resps)){

  # print(ParasiteTypes[q])
  #
  # r <- 1
  #
  # for(r in r:length(ParasiteTypes)){

  FocalResp <- Resps[q]

  print(FocalResp)

  iucn %>%
    select(FocalResp, LogRichness,
           X, Y,
           FullCovar, Sp) %>%
    na.omit %>% droplevels ->
    TestDF

  TestDF[,"Response"] <- TestDF[,FocalResp]

  TestDF %>%
    group_by(hOrder) %>%
    summarise(N = n(),
              Prevalence = Prev(Response),
              NInf = sum(Response>0)) %>%
    filter(N>20, Prevalence > 0.001, NInf > 1) %>%
    pull(hOrder) ->

    OrderInclude

  TestDF %>%
    filter(hOrder %in% OrderInclude) %>%
    droplevels ->

    TestDF

  print(nrow(TestDF))

  TestDF %>%
    select(FullCovar) %>%
    mutate_if(is.numeric, ~.x %>% scale %>% c) ->
    TestDF[,FullCovar]

  if(q < 2){

    IM1 <- INLAModelAdd(
      Response = FocalResp,
      Explanatory = FullCovar, # %>% setdiff("LogCites"), # %>% setdiff("UrbanBinary"),
      # Add = c(#"UrbanBinary:hOrder",
      #         # "LogCites",
      #         # "Traded:hOrder",
      #         #"UrbanBinary:LogCites"
      #   # "UrbanBinary", "UrbRurPopRatioLn"
      #   ),
      # Add = c("LogCites"),
      Family = FamilyList[[Resps[q]]],
      Data = TestDF,
      Base = T, Beep = F,
      AllModels = T,
      AddSpatial = T

    )

  }else{

    IM1 <- INLAModelAdd(
      Response = FocalResp,
      Explanatory = FullCovar, # %>% setdiff("LogCites"), # %>% setdiff("UrbanBinary"),
      # Add = c(#"UrbanBinary:hOrder",
      #         # "LogCites",
      #         # "Traded:hOrder",
      #         #"UrbanBinary:LogCites"
      #   # "UrbanBinary", "UrbRurPopRatioLn"
      #   ),
      Add = c("LogRichness"),
      Family = FamilyList[[Resps[q]]],
      Data = TestDF,
      Base = T, Beep = F,
      AllModels = T,
      AddSpatial = T

    )


  }


  # IM1$FinalModel %>% Efxplot

  IMList[[FocalResp]] <- IM1

  # IM1$FinalModel %>% Efxplot %>% plot

  # IMList[[FocalResp]] %>%
  #   saveRDS(file = paste0("Output/Models_", FocalResp, ".rds"))

}

library(patchwork)

IMList %>% map("FinalModel") %>% Efxplot + ggtitle("Base") +
  IMList %>% map(c("Spatial", "Model")) %>% Efxplot + ggtitle("Spatial") +
  plot_layout(guides = "collect")
