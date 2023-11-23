library(tidyverse)
library(ggExtra)
#library(matlib)
library(MASS)
library(purrr)

theme_set(theme_bw())

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

DBS_spectra <- read_tsv('/Volumes/Nexus/Users/busha/Analysis/2021-07-07-FOOOF-reparam/data_30kHz/ITI_power_spectra_30kHz.tsv') %>%
  mutate(electrode=tolower(electrode))

DBS_per <- read_csv('data/DBS_periodic_components_all.csv')%>%
  mutate(electrode=tolower(electrode)) 

#calculating gaussian amplitudes from peak heights
calc_gauss_amp <- function(p,...){
  rank <- sum(!is.na(p$PH1),!is.na(p$PH2),!is.na(p$PH3),!is.na(p$PH4),!is.na(p$PH5),!is.na(p$PH6))
  print(rank)
  
  A <- array(0,c(6,1))
  if (rank>0){
    M<-array(0,c(rank,rank))
    PW<-array(0,c(rank,1))
    for(i in 1:rank){
      for(j in 1:rank){
        M[i,j] <- as.numeric(exp(-((p[1,paste0('CF',i)]-p[1,paste0('CF',j)])^2/(2*(p[1,paste0('BW',j)]/2)^2))))
        PW[i] <- as.numeric(p[1,paste0('PH',i)])
      }
    }
    if(rank > 1){
      A[1:rank] <- ginv(M) %*% PW
    }else if(rank==1){
      A[1] <- PW[1]
    }
  }
  A <- as.data.frame(t(A))
  names(A) <- paste0('Gamp',1:6)
  return(as.data.frame(A))
}

DBS_per_Gamp <- DBS_per %>%
  group_by(subject,session_id,electrode) %>%
  group_modify(calc_gauss_amp,.keep=TRUE) 

DBS_per <- DBS_per %>%
  left_join(DBS_per_Gamp)

fmin = 1

STN_DBS_spectra <- DBS_spectra %>%
  filter((subject=='DBS3028' & session_id==4 & electrode %in% c('dbs_l1','ecog_248')) | #STN
         (subject=='DBS4067' & session_id==2 & electrode %in% c('dbs_r3','ecog_203')) | #VIM
         (subject=='DBS4079' & session_id==4 & electrode %in% c('dbs_l3a')) | #GPi
         (subject=='DBS4079' & session_id==3& electrode %in% c('ecog_153'))) %>% #GPi 
  pivot_longer(cols=6:last_col(),names_to = "frequency", values_to="power") %>%
  mutate(frequency=as.numeric(frequency)) %>%
  left_join(DBS_aper) %>%
  left_join(DBS_per) %>%
  mutate(comp_aper = offset + log10(10^(log_knee*exponent)+fmin^exponent) - log10(10^(log_knee*exponent)+frequency^exponent)) %>%
  mutate(comp_per1 = Gamp1 * exp(-((frequency-CF1)^2/(2*(BW1/2)^2)))) %>%
  mutate(comp_per2 = Gamp2 * exp(-((frequency-CF2)^2/(2*(BW2/2)^2)))) %>%
  mutate(comp_per3 = Gamp3 * exp(-((frequency-CF3)^2/(2*(BW3/2)^2)))) %>%
  mutate(comp_per4 = Gamp4 * exp(-((frequency-CF4)^2/(2*(BW4/2)^2)))) %>%
  mutate(comp_per5 = Gamp5 * exp(-((frequency-CF5)^2/(2*(BW5/2)^2)))) %>%
  mutate(comp_per6 = Gamp6 * exp(-((frequency-CF6)^2/(2*(BW6/2)^2)))) %>%
  replace_na(list(comp_per1=0,comp_per2=0,comp_per3=0,comp_per4=0,comp_per5=0,comp_per6=0)) %>%
  mutate(fooof = comp_aper + comp_per1 + comp_per2 ++ comp_per3 + comp_per4 + comp_per5 + comp_per6) %>%
  left_join(subjects) %>%
  mutate(electrode_type_dbs_target = paste(electrode_type,dbs_target,sep='_')) %>%
  mutate(electrode_type_dbs_target = factor(electrode_type_dbs_target,levels=c('dbs_VIM','ecog_VIM','dbs_STN','ecog_STN','dbs_GPi','ecog_GPi'))) %>%
  mutate(dbs_target=factor(dbs_target,levels=c('VIM','STN','GPi')))


color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG
color_ep_seeg = '#8A4F80' #EP sEEG
color_pd_stn = '#F7924A' #PD STN
color_pd_gpi = '#F9BD00' #PD GPi
color_et_vim = '#36A5D1' #ET VIM
color_ep_cm = '#9EB859' #EP CM
palette = c(color_et_vim,color_et_ecog,color_pd_stn,color_pd_ecog,color_pd_gpi,color_pd_ecog)


STN_DBS_spectra %>%
  mutate(dbs_target=factor(dbs_target,levels=c("STN","VIM","GPi"))) %>%
  filter(frequency>= fmin & frequency<=15000) %>%
  ggplot() +
  geom_line(aes(x=frequency,y=power,color=electrode_type_dbs_target),linetype=1)+
  geom_line(aes(x=frequency,y=comp_aper,group=electrode_type),
            data=STN_DBS_spectra %>% filter(frequency>= fmin & frequency<=250),
            linetype=1,color='gray')+
  geom_line(aes(x=frequency,y=fooof,group=electrode_type),
            data=STN_DBS_spectra %>% filter(frequency>= fmin & frequency<=250),          
            linetype=3,color='black')+
  scale_color_manual(values=palette) +
  facet_wrap(~dbs_target,ncol=1,scales='free')+
  scale_x_log10() +
  labs(x='Frequency (Hz)',y='Spectral Power Density')
ggsave('fig/B18_01_PSD-and-fits.pdf',width=8,height=9)



