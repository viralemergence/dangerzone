
# 03_Greg Analysis ####

rm(list = ls())

{

  library(tidyverse); library(data.table); library(magrittr); library(ggregplot); library(cowplot)
  library(geiger);library(ape);library(picante); library(patchwork); library(colorspace)
  library(tidyverse); library(INLA); library(ggregplot); library(glue); library(fs)

  theme_set(theme_cowplot())

}


###########################################################
# Data reconciliation between Kayla and Greg
urban <- readRDS("~/Github/dangerzone/data/IUCNDF.rds")
iucn <- read_csv("~/Github/dangerzone/Data/FixedPathMar72022.csv")
iucn %<>%
  mutate(Endangered = as.logical(Endangered),
         Decreasing = as.logical(Decreasing),
         DataDeficient = as.logical(DataDeficient)) %>%
  select(-X8) %>%
  mutate(NVirion = round(10^NVirion - 1),
         NZoon = round(10^NZoon - 1),
         Citations = round(10^Citations - 1))
urban %>%
  select(-c(Citations,NumVirus,NVirion,NZoon,Endangered,Decreasing,Data_Deficient)) %>%
  left_join(iucn) %>%
  mutate(LogCites = log(Citations + 1)) ->
        iucn

iucn$NVirion[iucn$NVirion==0] <- NA # I THINK this harmonizes it
iucn$NVirion[iucn$NZoon==0] <- NA
###########################################################

iucn %<>% filter(!is.na(DataDeficient))

# iucn %<>% rename(DataDeficient = Data_Deficient)

# source("0_Urban Import.R")

dir_create(c("Output", "Figures"))

# UrbanDF %<>% filter(!Pseudoabsence) # Greg filters

Resps <- c(

  "NVirion", "NZoon"

)

FamilyList <- rep("nbinomial", 2)

names(FamilyList) <- Resps

FullCovar <- c("hOrder",
               "LogCites",
               "Endangered", "DataDeficient", "Decreasing",
               # "populationTrend",
               "DomesticBinary",
               # "HumanDistance",
                "LogArea",
               # "DietDiversity",
               # "LogMass",
               "PC1", "PC2")

PhyloList <-
  IMList <-
  TestDFList <-
  list()

if(1){ # Removing pseudoabsences

  iucn %<>%
    filter(!is.na(NVirion)) %>% # Removing all pseudoabsences
    # filter(!(is.na(NVirion) & LogCites == 0)) %>% # Removing pseudoabsences and zero-citations
    # mutate_at("NZoon", ~ifelse(is.na(.x), 0, .x)) %>%
    # mutate_at("NVirion", ~.x - 1) %>%
    mutate_at(c("NZoon", "NVirion"), ~ifelse(is.na(.x), 0, .x))

  iucn %>% nrow %>% print

}else{ # Including pseudoabsences

  iucn %<>%
    mutate_at(c("NZoon", "NVirion"), ~ifelse(is.na(.x), 0, .x))

  iucn %>% nrow %>% print

}

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

  # TestDF %<>% mutate_at("populationTrend", ~factor(.x, levels = c("Stable", "Decreasing", "Increasing", "Unknown")))

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
      AddSpatial = T,
      Groups = T, GroupVar = "Endangered"

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
      AddSpatial = T,
      Groups = T, GroupVar = "Endangered"

    )


  }


  # IM1$FinalModel %>% Efxplot

  IMList[[FocalResp]] <- IM1

  # IM1$FinalModel %>% Efxplot %>% plot

  IMList[[FocalResp]] %>%
    saveRDS(file = paste0("Output/Models_", FocalResp, ".rds"))

}

IMList %>% map("FinalModel") %>% Efxplot + ggtitle("Base") +
  IMList %>% map(c("Spatial", "Model")) %>% Efxplot + ggtitle("Spatial") +
  plot_layout(guides = "collect")

IMList %>% map("FinalModel") %>% Efxplot + ggtitle("Base") +
  IMList %>% map(c("Spatial", "Model")) %>% Efxplot + ggtitle("Spatial") +
  IMList %>% map(c("Spatial", "SpatiotemporalModel")) %>% Efxplot + ggtitle("Spatial") +
  plot_layout(guides = "collect")

# Figure ####

ModelEffects <-
  IMList %>%
  map(c("Spatial", "Model")) %>%
<<<<<<< HEAD
  Efxplot +
  # theme(legend.position = "top") +
=======
  Efxplot(PointSize = 3) +
  theme(legend.position = c(0.7, 0.07)) +
>>>>>>> df647dc12e34e94bb68b364d40abeb1ee551636e
  # theme(legend.position = "none") +
  scale_colour_manual(values = c(AlberColours[[2]], AlberColours[[3]]),
                      labels = c("Total",
                                 "Zoonotic"))



ModelEffects +  # Always check this
  scale_x_discrete(labels = rev(c("Intercept",
                              "Carnivora",
                              "Chiroptera",
                              "Primates",
                              "Rodentia",
                              "Citations",
                              "Endangered",
                              "Data Deficient",
                              "Domesticated",
                              "Geographic Range Area",
                              "Life History PC1",
                              "Life History PC2",
                              "Total Richness")))




Maps <-
  IMList %>%
  map(function(a){

    ggField(Model = a$Spatial$Model,
            Mesh = a$Spatial$Mesh, #Groups = 2,
            PointAlpha = 0.2,
            Points = a$Data[,c("X", "Y")]) +
      scale_fill_discrete_sequential(palette = AlberPalettes[[2]]) +
      theme_void()

  }) %>%
  ArrangeCowplot() +
  plot_layout(nrow = 2,
              guides = "collect")

Maps[[1]] <-
  Maps[[1]] + labs(fill = "Rich.")

Maps[[2]] <-
  Maps[[2]] + labs(fill = "Zoo.Rich.")

Maps[[2]] <- Maps[[2]] + scale_fill_discrete_sequential(palette = AlberPalettes[[3]])

(ModelEffects/Maps)# + plot_layout(widths = c(1, 6))

ModelEffects +
  theme(legend.position = c(0.7, 0.3))

ggsave("Figures/ModelPanels.jpeg", units = "mm", height = 200, width = 200)

# Path Analysis ####

(MeanEndangeredNV <- IMList$NVirion$Spatial$Model %>% GetEstimates("EndangeredTRUE"))
(MeanEndangeredNZ <- IMList$NZoon$Spatial$Model %>% GetEstimates("EndangeredTRUE"))
(MeanNVNZ <- IMList$NZoon$Spatial$Model %>% GetEstimates("LogRichness"))

(PEndangeredNV <- IMList$NVirion$Spatial$Model %>% INLAPValue("EndangeredTRUE") %>% unlist)
(PEndangeredNZ <- IMList$NZoon$Spatial$Model %>% INLAPValue("EndangeredTRUE") %>% unlist)
(PNVNZ <- IMList$NZoon$Spatial$Model %>% INLAPValue("LogRichness") %>% unlist)

Path <-
  data.frame(

    Var = c("Endangered_NVirus", "Endangered_NZoo", "NVirus_NZoo"),
    Estimates = c(MeanEndangeredNV, MeanEndangeredNZ, MeanNVNZ),
    P = c(PEndangeredNV, PEndangeredNZ, PNVNZ)

  )

Path %>% write.csv("Figures/PathCoefficients.csv", row.names = F)

EndangeredNV <- IMList$NVirion$Spatial$Model %>% GetEstimates("EndangeredTRUE", NDraws = 1000, Draw = T)
EndangeredNZ <- IMList$NZoon$Spatial$Model %>% GetEstimates("EndangeredTRUE", NDraws = 1000, Draw = T)
NVNZ <- IMList$NZoon$Spatial$Model %>% GetEstimates("LogRichness", NDraws = 1000, Draw = T)

Indirect <- EndangeredNV*NVNZ

Indirect %>% qplot

IndirectCI <- Indirect %>% as.mcmc %>% HPDinterval()
IndirectMean <- Indirect %>% as.mcmc %>% mean()

Indirect %>% qplot + geom_vline(aes(xintercept = c(IndirectCI, IndirectMean)), colour = "red")

# Endangered Maps ####

IMList$NZoon$Spatial$SpatiotemporalModel %>%
  ggField(Mesh = IMList$NZoon$Spatial$Mesh, Groups = 2,
          GroupVar = "Endangered",
          PointSize = 1,
          Points = IMList$NZoon$Data[,c("X", "Y", "Endangered")]) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]])


IMList$NZoon$Spatial$Model %>%
  ggField(Mesh = IMList$NZoon$Spatial$Mesh, #Groups = 2,
          Points = IMList$NZoon$Data[,c("X", "Y")]) +
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]])




