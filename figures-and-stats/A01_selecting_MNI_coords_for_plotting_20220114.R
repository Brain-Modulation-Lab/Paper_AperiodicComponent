library(tidyverse)

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

elec_session <- read_tsv('data/electrode_session.txt') %>%
  dplyr::select(subject,session_id,electrode,type,side,
         mni_nonlinear_x,mni_nonlinear_y,mni_nonlinear_z,target,
         fs_anatomy,DISTAL_label_1,DISTAL_label_2,DISTAL_label_3,MOREL_label_1,MOREL_label_2,MOREL_label_3,HCPMMP1_label_1) %>%
  mutate(electrode=tolower(electrode)) 


subjects<-read_tsv('data/subjects.tsv')

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
  #dplyr::select(-X1) %>%
  mutate(electrode=tolower(electrode)) %>%
  left_join(elec_session) %>%
  #filter(!is.na(mni_nonlinear_x)) %>%
  left_join(subjects %>% dplyr::select(subject,dx,dbs_target)) %>%
  filter(tolower(type) %in% c('ecog','dbs')) %>%
  mutate(target=tolower(target))

#removing contacts far from target structure
#overall 9.4% of dbs contacts are removed
DBS_aper_filt <- DBS_aper %>%
  rowwise() %>%
  mutate(on_gpi=any(c(DISTAL_label_1,DISTAL_label_2,DISTAL_label_3) %in%
                                     paste0(c('GPi_sensorimotor','GPi_postparietal','GPi_premotor'),rep(c('_L','_R'),each=3)))) %>%
  mutate(on_ventral_thalamus=any(c(MOREL_label_1,MOREL_label_2) %in%
                            paste0(c('VPM','VM','VLpv','VLa','VPLa','VPLp','VPLa','VPI','VLpd','VApc'),rep(c('_L','_R'),each=10)))) %>%
  mutate(on_stn=any(c(MOREL_label_1,MOREL_label_2) %in% c('STh_L','STh_R'))) %>%
  ungroup() %>%
  filter(!(type=='dbs' & target=='stn' & !on_stn)) %>%
  filter(!(type=='dbs' & target=='vim' & !on_ventral_thalamus)) %>%
#  filter(!(type=='dbs' & target=='gpi' & !on_gpi)) %>%
  dplyr::select(-on_gpi,-on_ventral_thalamus,-on_stn)

write_tsv(DBS_aper_filt,'data/A01_DBS_aper_coord_dx.tsv')

