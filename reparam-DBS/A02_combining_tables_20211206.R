library(tidyverse)
library(glue)
library(readr)
library(magrittr)

#PATH_ANALYSIS = "Z:/Users/busha/Analysis/2021-07-07-FOOOF-re[aram"
PATH_ANALYSIS = "/Volumes/Nexus/Users/busha/Analysis/2021-07-07-FOOOF-reparam"
PATH_DB = paste0(PATH_ANALYSIS,"/data")

setwd(PATH_ANALYSIS)

# Loading fooof results
subjects<-
  tibble(subject=paste0('DBS',c(seq(3001,3032),4057,4058,4061,seq(4063,4087)))) %>%
  mutate(path_spectra=glue(paste0(PATH_DB,'/{subject}_power_spectra.tsv'))) %>%
  mutate(exists_spectra=file.exists(path_spectra)) %>%
  filter(exists_spectra)

spectra <- subjects$path_spectra%>% 
  map_df(read_tsv)

spectra_w <- spectra %>%
  mutate(power=log10(power)) %>%
  filter(type=='iti') %>%
  filter(!str_detect(electrode,"^macro"))%>%
  filter(!str_detect(electrode,"^micro"))%>%
  group_by(subject,session_id) %>%
  mutate(is_dbs_session = any(str_detect(electrode,"^dbs"))) %>%
  ungroup() %>%
  #filter(is_dbs_session) %>%
  #select(-is_dbs_session) %>%
  pivot_wider(id_cols=c(subject,session_id,electrode,type,is_dbs_session),values_from=power,names_from=frequency)

spectra_w %>%
  write_tsv(glue('{PATH_DB}/ITI_power_spectra.tsv'))
