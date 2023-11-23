library(MASS)
library(tidyverse)
library(multidplyr)
library(ggExtra)
library(purrr)

theme_set(theme_bw())

clusters <- new_cluster(parallel::detectCores() - 2)

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')

artifacts <- read_tsv('data/artifacts.txt') %>%
  fill(subject,session_id,type) %>%
  mutate(subject=str_replace(subject,fixed('-'),'')) %>%
  mutate(electrode=tolower(electrode)) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) 

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
  select(-`...1`) %>%
  mutate(electrode=tolower(electrode)) 

DBS_spectra <- read_tsv('/Volumes/Nexus/Users/busha/Analysis/2021-07-07-FOOOF-reparam/data/ITI_power_spectra.tsv') %>%
  mutate(electrode=tolower(electrode))

DBS_per <- read_csv('data/DBS_periodic_components_all.csv')%>%
  mutate(electrode=tolower(electrode)) 


DBS_aper_coord_dx <- read_tsv('data/A01_DBS_aper_coord_dx.tsv') %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel)

MMP1_regions <- read_tsv('data/MMP1-regions2.tsv')

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG


##### ECoG

DBS_ecog_aper <- DBS_aper_coord_dx %>%
  left_join(MMP1_regions) %>%
  filter(electrode_type=='ecog') %>%
  filter(rsquared > 0.98) %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(!is.na(dbs_target)) %>%
  mutate(dx_dbs_target=str_c(dx,'_',dbs_target)) %>%
  mutate(dx_dbs_target=factor(dx_dbs_target,levels=c('PD_STN','PD_GPi','ET_VIM'),ordered=TRUE)) %>%
  mutate(sub_ses_el = str_c(subject,'_',session_id,'_',electrode))

DBS_ecog_spectra <- DBS_spectra  %>%
  mutate(sub_ses_el = str_c(subject,'_',session_id,'_',electrode),.after=electrode)%>%
  filter(sub_ses_el %in% DBS_ecog_aper$sub_ses_el) %>%
  pivot_longer(cols=7:last_col(),names_to = "frequency", values_to="power") %>%
  mutate(frequency=as.numeric(frequency)) %>%
  left_join(DBS_ecog_aper %>% select(subject, session_id, electrode, region1,HCPMMP1_label_1,dx,dbs_target,dx_dbs_target))

DBS_ecog_spectra_agg_region1 <- DBS_ecog_spectra %>%
  group_by(region1,dx_dbs_target,frequency) %>%
  summarise(mean_power = mean(power, na.rm=TRUE), 
            sd_power = sd(power, na.rm=TRUE), 
            mad_power = mad(power, na.rm=TRUE), 
            median_power = median(power, na.rm=TRUE), 
            n=n()) %>%
  mutate(sem=sd_power/sqrt(n))
  
DBS_ecog_spectra_agg_region1 %>%
  ggplot() +
  aes(y=median_power,x=frequency,ymax=median_power+mad_power,ymin=median_power-mad_power,color=dx_dbs_target) +
  geom_line() +
  geom_smooth(stat="identity") +
  scale_x_log10() +
  facet_wrap(~region1)


color_pd_stn_ecog = '#6F67A6' 
color_pd_gpi_seeg = '#8A4F80' 
color_et_ecog = '#C4604F' 
color_pd_stn = '#F7924A' #PD STN
color_pd_gpi = '#F9BD00' #PD GPi
color_et_vim = '#36A5D1' #ET VIM

DBS_ecog_spectra_agg_MMP1 <- DBS_ecog_spectra %>%
  mutate(sub_el = str_c(subject,'_',electrode)) %>%
  group_by(region1,HCPMMP1_label_1,dx,dbs_target,dx_dbs_target,frequency) %>%
  #partition(clusters) %>%
  summarise(mean_power = mean(power, na.rm=TRUE), 
            sd_power = sd(power, na.rm=TRUE), 
            mad_power = mad(power, na.rm=TRUE), 
            median_power = median(power, na.rm=TRUE), 
            n=n(),
            n_sub = length(unique(subject)),
            n_el = length(unique(sub_el))) %>%
  mutate(sem=sd_power/sqrt(n)) %>%
  mutate(region_MMP1 = str_c(region1,' [',HCPMMP1_label_1,']')) %>%
  filter(!is.na(region1)) %>%
  mutate(label = str_c(dx,'[',dbs_target,'] Ns=',n_sub,' Ne=',n_el)) %>%
  write_tsv("fig/B18b_01_ECoG-group-level-MMP1.tsv")

DBS_ecog_spectra_agg_MMP1 %>%
  ggplot() +
  aes(y=median_power,x=frequency,ymax=median_power+mad_power,ymin=median_power-mad_power,color=dx_dbs_target) +
  geom_line() +
  geom_smooth(stat="identity") +
  geom_text(aes(label=label,y=-1*as.numeric(dx_dbs_target),x=NULL,fill=NULL),x=0,hjust=0,parse=FALSE,size=3) +
  scale_y_continuous(breaks=-3:3,labels=10^(-3:3)) +
  scale_x_log10() +
  facet_wrap(~region_MMP1,ncol=6) +
  scale_color_manual(values=c(color_pd_stn_ecog,color_pd_gpi_seeg,color_et_ecog)) +
  theme(legend.position='none') +
  coord_cartesian(ylim=c(-4,4),xlim=c(1,250))+
  labs(x='Frequency (Hz)', y="Power (uV^2/Hz)")
 ggsave("fig/B18b_01_ECoG-group-level-MMP1.pdf",width=16,height=18)

 
 
 ##### Subcortical
 
 artifacts <- read_tsv('data/artifacts.txt') %>%
   fill(subject,session_id,type) %>%
   mutate(subject=str_replace(subject,fixed('-'),'')) %>%
   mutate(electrode=tolower(electrode)) %>%
   mutate(subsesel = paste(subject,session_id,electrode,sep='_')) 
 
 DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
   select(-'...1') %>%
   mutate(electrode=tolower(electrode)) %>%
   left_join(subjects %>% select(subject,dx,dbs_target)) %>%
   mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
   filter(!subsesel %in% artifacts$subsesel) %>%
   filter(!str_detect(electrode,fixed('dbs_r'))) 
 
 DBS_aper_coord_dx <- read_tsv('data/A01_DBS_aper_coord_dx.tsv') %>%
   mutate(electrode=tolower(electrode)) %>%
   mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
   filter(!subsesel %in% artifacts$subsesel) %>%
   filter(!str_detect(electrode,fixed('dbs_r'))) 
 
 DBS_LFP_aper <- DBS_aper_coord_dx %>%
   filter(electrode_type=='dbs') %>%
   filter(dx%in%c('PD','ET')) %>%
   filter(!is.na(dbs_target)) %>%
   mutate(dx_dbs_target=str_c(dx,'_',dbs_target)) %>%
   mutate(dx_dbs_target=factor(dx_dbs_target,levels=c('PD_STN','PD_GPi','ET_VIM'),ordered=TRUE)) %>%
   mutate(sub_ses_el = str_c(subject,'_',session_id,'_',electrode))
 
 DBS_LFP_spectra <- DBS_spectra  %>%
   mutate(sub_ses_el = str_c(subject,'_',session_id,'_',electrode),.after=electrode)%>%
   filter(sub_ses_el %in% DBS_LFP_aper$sub_ses_el) %>%
   pivot_longer(cols=7:last_col(),names_to = "frequency", values_to="power") %>%
   mutate(frequency=as.numeric(frequency)) %>%
   left_join(DBS_LFP_aper %>% select(subject, session_id, electrode, dx,dbs_target,dx_dbs_target))
 
 DBS_LFP_spectra_agg_target <- DBS_LFP_spectra %>%
   mutate(sub_el = str_c(subject,'_',electrode)) %>%
   group_by(dx,dbs_target,dx_dbs_target,frequency) %>%
   summarise(mean_power = mean(power, na.rm=TRUE), 
             sd_power = sd(power, na.rm=TRUE), 
             mad_power = mad(power, na.rm=TRUE), 
             median_power = median(power, na.rm=TRUE), 
             n=n(),
             n_sub = length(unique(subject)),
             n_el = length(unique(sub_el))) %>%
   mutate(sem=sd_power/sqrt(n)) %>%
   mutate(label = str_c(dx,'[',dbs_target,'] Ns=',n_sub,' Ne=',n_el)) %>%
   write_tsv("fig/B18b_02_LFP-group-level-DBS-target.tsv")

 color_pd_stn = '#F7924A' #PD STN
 color_pd_gpi = '#F9BD00' #PD GPi
 color_et_vim = '#36A5D1' #ET VIM
 
 DBS_LFP_spectra_agg_target %>%
   ggplot() +
   aes(y=median_power,x=frequency,ymax=median_power+mad_power,ymin=median_power-mad_power,color=dx_dbs_target) +
   geom_line() +
   geom_smooth(stat="identity") +
   geom_text(aes(label=label,x=NULL),x=0,y=-2,hjust=0,parse=FALSE,size=3) +
   scale_y_continuous(breaks=-3:3,labels=10^(-3:3)) +
   scale_x_log10() +
   facet_wrap(~dbs_target,ncol=3) +
   scale_color_manual(values=c(color_pd_gpi,color_pd_stn,color_et_vim)) +
   theme(legend.position='none') +
   coord_cartesian(ylim=c(-4,4),xlim=c(1,250))+
   labs(x='Frequency (Hz)', y="Power (uV^2/Hz)")
 ggsave("fig/B18b_02_LFP-group-level-DBS-target.pdf",width=8.2,height=2.5)
 


