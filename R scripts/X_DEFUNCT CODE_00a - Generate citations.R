library(easyPubMed)
library(magrittr)
library(tidyverse)

CiteCounter <- function(Name) {
  Name %>% str_trim %>%
    get_pubmed_ids %>%
    extract2("Count") %>%
    as.character %>%
    as.numeric
}

iucn %<>% select(Host)
iucn$Citations = NA
i <- 1

for(i in i:nrow(iucn)){
  Sp <- iucn$Host[i]
  print(Sp)
  iucn$Citations[i] <- CiteCounter(Sp)
  print(paste0("Citations: ", iucn$Citations[i]))
}

write_csv(iucn, "Citation2.csv")