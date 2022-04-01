
# 0_Urban Import #####

{

  library(tidyverse); library(data.table); library(magrittr); library(ggregplot); library(cowplot)
  library(geiger);library(ape);library(picante)

  theme_set(theme_cowplot())

}

iucn %<>% mutate(Sp = Host %>% str_replace_all(" ", "_") %>% CamelConvert)

Panth1 <- read.delim("UOData/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename_all(~str_replace(.x, "MSW05_", "h")) %>%
  rename(Sp = hBinomial) %>%
  mutate_at("Sp", ~.x %>% str_trim %>% str_replace_all(" ", "_"))

iucn %<>% left_join(Panth1)

Panth1 %>%
  filter(hOrder%in%c("Cetacea", "Sirenia")|
           hFamily%in%c("Phocidae", "Odobenidae", "Otariidae")) %>%
  pull(Sp) ->
  MarineSp

Domestic <- read.csv("UOData/Domestic.csv")

Domestic$Species.and.subspecies %>%
  str_split("[(]") -> Bracket1 # Detecting species names

Bracket1 %>% map(~.x[str_detect(.x, "[)]")]) -> Bracket2 # detecting incomplete species names

Bracket2 %>% map(~.x %>% str_split("[)]") %>% map(1)) -> Bracket3

1:length(Bracket3) %>% map(function(b){

  # print(b)

  a <- Bracket3[[b]]

  Genus <- a[[1]] %>% str_split(" ") %>% extract2(1) %>% extract2(1)

  if((a[[1]] %>% str_split(" ") %>% extract2(1) %>% length)>1){

    a %>% map(~.x %>% str_split(" ") %>% map_chr(2) %>% paste0(Genus, "_", .))

  }else NA

}) %>% unlist %>% unique %>% sort -> DomesticSpecies

AddDomestic <- c("Capra_hircus",
                 "Neovison_vison")

DomesticSpecies %<>% c(AddDomestic) %>% unique %>% sort

iucn %>%
  mutate(DomesticBinary = as.numeric(Sp %in% DomesticSpecies) %>% as.factor) ->
  iucn

CitationDF <- readRDS("Data/CitationDF.rds")

iucn %<>% left_join(CitationDF, by = "Sp") %>%
  mutate(LogCites = log10(Citations + 1))

# Phylogenetic data ####

STFull <- read.nexus("UOData/ele_1307_sm_sa1.tre")[[1]]
FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix

FullSTMatrix["Homo_sapiens",] %>% as.data.frame() %>%
  rownames_to_column %>% rename("HumanDistance" = ".") %>%
  left_join(iucn, ., by = c("Sp" = "rowname")) ->
  iucn

# Fixing phenotypic stuff ####

iucn %<>%
  rename(Mass = X5.1_AdultBodyMass_g,
         HomeRange = X22.1_HomeRange_km2) %>%
  mutate(LogMass = log(Mass))

DegreeGet <- function(a){

  a[a>0] <- 1
  rowSums(a) %>% return

}

# Adding Spatial Data ####

load("UOData/FullRangeOverlapMercator.Rdata") # From Albersnet ####

FullRangeAdj %>% colSums -> SympatryStrengthVector
FullRangeAdj %>% DegreeGet -> SympatryDegreeVector

SympatryStrengthVector[iucn$Sp] -> iucn$SympatryStrength
SympatryDegreeVector[iucn$Sp] -> iucn$SympatryDegree

iucn %>%
  select(contains("Sympatry")) %>%
  mutate_all(log10) ->
  iucn[,paste0("Log", c("SympatryStrength", "SympatryDegree"))]

# Importing Eskew Data ####

EskewPCA <- read.csv("EskewFiles/reservoir_data_for_greg.csv")

EskewPCA %<>% mutate_at("PC1", ~ -.x)

EskewPCA %>%
  mutate(Sp = binom %>% str_trim %>%  str_replace_all(" ", "_")) %>%
  select(PC1:Sp) %>%
  left_join(iucn, ., by = "Sp") ->
  iucn

# Eliminating non-Eutherians ####

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia",
                   "Notoryctemorphia",
                   "Monotremata")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

iucn %<>% filter(!Sp %in% NonEutherianSp)

iucn %>% dim

# Adding centroids etc ####

CentroidList <- readRDS("UOData/CentroidList.rds")
ContinentsInhabitedList <- readRDS("UOData/ContinentsInhabitedList.rds")

ContinentsInhabitedList %>% unlist %>% unique -> ContinentVars

CentroidList %<>%
  map(~.x %>%
        matrix(nrow = 2) %>% t %>% as.data.frame %>% rename(X = V1, Y = V2) %>%
        summarise_all(mean)
  )

ContinentDF <-
  data.frame(Sp = names(ContinentsInhabitedList))

ContinentDF$Continents <- ContinentsInhabitedList

ContinentDF %<>%
  unnest(Continents) %>%
  mutate(Presence = 1) %>%
  tidyr::pivot_wider(#id_cols = "Sp",
    names_from = "Continents", values_from = "Presence") %>%
  mutate_at(-1, ~ifelse(is.na(.x), 0, .x))

CentroidList %>% bind_rows(.id = "Sp") %>%
  left_join(ContinentDF, .) %>%
  left_join(iucn, .) -> iucn

iucn %>%
  select(Eurasia:`NAm`) %>% names -> ContinentVar

iucn %>%
  dplyr::select(all_of(ContinentVar)) %>% rowSums ->
  iucn[,"NContinents"]

iucn %<>% mutate_at("NContinents", ~ifelse(.x == 0, 1, .x))

AreaList <- read.csv("UOData/RangeAreas.csv")

iucn %<>%

  left_join(AreaList, by = c("Sp" = "Species"))

(iucn$Area+1) %>%
  log -> iucn$LogArea

# Elton ####

Elton <- read.delim("UOData/Defunct/MamFuncDat.txt")

Elton %<>%
  mutate(Sp = Scientific %>% str_trim %>% str_replace_all(" ", "_"))

Elton %>%
  dplyr::select(contains("Diet.")) %>%
  mutate_all(AsBinary) %>%
  rowSums(na.rm = T) -> Elton$DietNo

Elton %>% dplyr::select(Sp, DietNo) %>% left_join(iucn, .) -> iucn

library(vegan)

Elton %>%
  dplyr::select(contains("Diet.")) %>%
  dplyr::select(-matches("Source|Certainty")) -> Proportions

Proportions %>% diversity -> Diversities

Elton$DietDiversity <- Diversities

# Elton %>% ggplot(aes(DietNo, DietDiversity)) + geom_point() + geom_smooth(method = lm)

Elton %>% dplyr::select(Sp, DietDiversity) %>%
  left_join(iucn, .) ->
  iucn

iucn %<>%
  mutate(Endangered = str_detect(redlistCategory, "Endangered")) %>%
  mutate(Data_Deficient = str_detect(redlistCategory, "Data"))

# Adding population trend ####

PopTrend <- read.csv("Data/PopulationTrend.csv")

iucn <-
  PopTrend %>% mutate(Sp = Host %>% str_replace(" ", "_")) %>%
  dplyr::select(Sp, populationTrend) %>%
  left_join(iucn, .)

iucn %<>%
  mutate(Decreasing = (populationTrend == "Decreasing"))

iucn %>% saveRDS("IUCNDF.rds")
