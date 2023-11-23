library(tidyverse)
library(ggExtra)
library(matlib)
library(purrr)

theme_set(theme_bw())

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

fmin = 1

elecLoc <- read_tsv('~/Dropbox (Personal)/Lab-BML/Expan/2021-04-01-CM-vs-cortex-elecLoc/data/electrode-localizations.tsv') %>%
  rename(subject=Subject) %>%
  mutate(subject=str_replace(subject,'EM','sEEG-')) %>%
  mutate(electrode=tolower(label)) %>%
  mutate(subel = paste(subject,electrode,'_'))
  
sEEG_aper <- read_csv('data/sEEG_aperiodic_components_all.csv') %>%
  select(-X1) %>%
  mutate(electrode=tolower(electrode)) %>%
  select(-fooof_export_path) 

sEEG_per <- read_csv('data/sEEG_periodic_components_all.csv')%>%
  select(-X1) %>%
  mutate(electrode=tolower(electrode)) %>%
  select(-fooof_export_path)

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
      A[1:rank] <- inv(M) %*% PW
    }else if(rank==1){
      A[1] <- PW[1]
    }
  }
  A <- as.data.frame(t(A))
  names(A) <- paste0('Gamp',1:6)
  return(as.data.frame(A))
}

sEEG_per <- read_csv('data/sEEG_periodic_components_all.csv')%>%
  select(-X1) %>%
  mutate(electrode=tolower(electrode)) %>%
  select(-fooof_export_path) 

sEEG_per_Gamp <- sEEG_per %>%
  group_by(subject,session_id,electrode,localization) %>%
  group_modify(calc_gauss_amp,.keep=TRUE) 
  
sEEG_per <- sEEG_per %>%
  left_join(sEEG_per_Gamp)


sEEG_spectra <- read_csv('data/sEEG_PSD_hp_fir_1hz.csv') %>%
  select(-X1) %>%
  mutate(electrode=tolower(electrode)) %>%
  pivot_longer(cols=6:394,names_to = "frequency", values_to="power") %>%
  mutate(frequency=as.numeric(frequency)) %>%
  mutate(subel = paste(subject,electrode,'_')) %>%
  left_join(sEEG_aper) %>%
  left_join(sEEG_per) %>%
  mutate(comp_aper = offset + log10(10^(log_knee*exponent)+fmin^exponent) - log10(10^(log_knee*exponent)+frequency^exponent)) %>%
  mutate(comp_per1 = Gamp1 * exp(-((frequency-CF1)^2/(2*(BW1/2)^2)))) %>%
  mutate(comp_per2 = Gamp2 * exp(-((frequency-CF2)^2/(2*(BW2/2)^2)))) %>%
  mutate(comp_per3 = Gamp3 * exp(-((frequency-CF3)^2/(2*(BW3/2)^2)))) %>%
  mutate(comp_per4 = Gamp4 * exp(-((frequency-CF4)^2/(2*(BW4/2)^2)))) %>%
  mutate(comp_per5 = Gamp5 * exp(-((frequency-CF5)^2/(2*(BW5/2)^2)))) %>%
  mutate(comp_per6 = Gamp6 * exp(-((frequency-CF6)^2/(2*(BW6/2)^2)))) %>%
  replace_na(list(comp_per1=0,comp_per2=0,comp_per3=0,comp_per4=0,comp_per5=0,comp_per6=0)) %>%
  mutate(fooof = comp_aper + comp_per1 + comp_per2 + comp_per3 + comp_per4 + comp_per5 + comp_per6) 

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG
color_ep_seeg = '#C4604F' #EP sEEG
color_pd_stn = '#F7924A' #PD STN
color_pd_gpi = '#F9BD00' #PD GPi
color_et_vim = '#36A5D1' #ET VIM
color_ep_cm = '#36A5D1' #EP CM
palette = c(color_ep_seeg,color_ep_cm)


sEEG_spectra %>%
  filter(subject=='sEEG-1023' & electrode %in% c('lcm3','lpi7') & session_id==0) %>%
  filter(frequency>= fmin & frequency<=250) %>% 
  ggplot() +
  geom_line(aes(x=frequency,y=power,color=localization),linetype=1)+
  geom_line(aes(x=frequency,y=comp_aper,group=localization),linetype=1,color='gray')+
  geom_line(aes(x=frequency,y=fooof,group=localization),linetype=3,color='black')+
  scale_color_manual(values=palette) +
  facet_wrap(~subject,ncol=1,scales='free')+
  scale_x_log10() +
  labs(x='Frequency (Hz)',y='Spectral Power Density')
ggsave('fig/B09_01_PSD-and-fits.pdf',width=5,height=3)

