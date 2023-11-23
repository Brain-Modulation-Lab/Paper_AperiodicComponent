library(tidyverse)
library(ggExtra)
library(coin)

theme_set(theme_bw())

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')

sEEG_aper <- read_csv('data/sEEG_aperiodic_components_all.csv') %>%
  select(-...1) %>%
  mutate(session_id=ceiling(session_id)) %>%
  mutate(electrode=tolower(electrode)) %>%
  mutate(electrode_type=factor(localization,levels=c("subcortical","cortical"))) %>%
  mutate(offset = offset + 12) #units are in uV, but interpreted as being in Volts

se <- function(x) sd(x)/sqrt(length(x))

sEEG_aper_agg2 <- sEEG_aper %>%
  group_by(subject,electrode_type) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
            sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
            offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) 

#====


artifacts <- read_tsv('data/artifacts.txt') %>%
  fill(subject,session_id,type) %>%
  mutate(subject=str_replace(subject,fixed('-'),'')) %>%
  mutate(electrode=tolower(electrode)) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) 

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
  select(-...1) %>%
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


se <- function(x) sd(x)/sqrt(length(x))

DBS_aper_agg2 <- DBS_aper_coord_dx %>%
  filter(tolower(type)%in%c('ecog','dbs'))%>%
  filter(toupper(dbs_target)%in%c('STN','GPI','VIM')) %>%
  filter(toupper(dx)%in%c('PD','ET')) %>%
  group_by(subject,dbs_target,electrode_type) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
            sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
            offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) %>%
  mutate(electrode_type_dbs_target = paste(electrode_type,dbs_target,sep='_')) %>%
  mutate(electrode_type_dbs_target = factor(electrode_type_dbs_target,levels=c('dbs_VIM','ecog_VIM','dbs_STN','ecog_STN','dbs_GPi','ecog_GPi'))) %>%
  mutate(dbs_target=factor(dbs_target,levels=c('VIM','STN','GPi')))

#====

aper_agg2 <- DBS_aper_agg2 %>% 
  ungroup() %>%
  mutate(cohort = str_c(dx,'_',dbs_target,'_DBS')) %>%
  mutate(is_cortical = electrode_type == "ecog") %>%
  select(subject,cohort,is_cortical,offset,log_knee,exponent) %>%
  bind_rows(
    sEEG_aper_agg2 %>% 
      ungroup() %>%
      mutate(cohort = 'EP_CM_sEEG') %>%
      mutate(is_cortical = electrode_type == "cortical") %>%
      select(subject,cohort,is_cortical,offset,log_knee,exponent))
write_tsv(aper_agg2,'data/aper-median-by-subject-all-cohorts.tsv')

#checking consistency
aper_agg2 %>%
group_by(cohort, is_cortical) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
            sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
            offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  View()



pw <- tibble()
cohorts = unique(aper_agg2$cohort)
for(i in 1:3){
  for (j in (i+1):4){
    pval <- pvalue(oneway_test(exponent ~ factor(cohort), 
                        data=aper_agg2 %>% 
                          filter(is_cortical) %>% 
                          filter(cohort %in% c(cohorts[i],cohorts[j])), distribution=approximate(nresample=9999)))
    pw <- bind_rows(pw,tibble(param='exponent',c1=cohorts[i],c2=cohorts[j],pvalue=as.numeric(pval)))
  }
}
for(i in 1:3){
  for (j in (i+1):4){
    pval <- pvalue(oneway_test(offset ~ factor(cohort), 
                               data=aper_agg2 %>% 
                                 filter(is_cortical) %>% 
                                 filter(cohort %in% c(cohorts[i],cohorts[j])), distribution=approximate(nresample=9999)))
    pw <- bind_rows(pw,tibble(param='offset',c1=cohorts[i],c2=cohorts[j],pvalue=as.numeric(pval)))
  }
}
for(i in 1:3){
  for (j in (i+1):4){
    pval <- pvalue(oneway_test(log_knee ~ factor(cohort), 
                               data=aper_agg2 %>% 
                                 filter(is_cortical) %>% 
                                 filter(cohort %in% c(cohorts[i],cohorts[j])), distribution=approximate(nresample=9999)))
    pw <- bind_rows(pw,tibble(param='log_knee',c1=cohorts[i],c2=cohorts[j],pvalue=as.numeric(pval)))
  }
}
pw_cortical <- pw %>% mutate(is_cortical=1)


pw <- tibble()
cohorts = unique(aper_agg2$cohort)
for(i in 1:3){
  for (j in (i+1):4){
    pval <- pvalue(oneway_test(exponent ~ factor(cohort), 
                               data=aper_agg2 %>% 
                                 filter(!is_cortical) %>% 
                                 filter(cohort %in% c(cohorts[i],cohorts[j])), distribution=approximate(nresample=9999)))
    pw <- bind_rows(pw,tibble(param='exponent',c1=cohorts[i],c2=cohorts[j],pvalue=as.numeric(pval)))
  }
}
for(i in 1:3){
  for (j in (i+1):4){
    pval <- pvalue(oneway_test(offset ~ factor(cohort), 
                               data=aper_agg2 %>% 
                                 filter(!is_cortical) %>% 
                                 filter(cohort %in% c(cohorts[i],cohorts[j])), distribution=approximate(nresample=9999)))
    pw <- bind_rows(pw,tibble(param='offset',c1=cohorts[i],c2=cohorts[j],pvalue=as.numeric(pval)))
  }
}
for(i in 1:3){
  for (j in (i+1):4){
    pval <- pvalue(oneway_test(log_knee ~ factor(cohort), 
                               data=aper_agg2 %>% 
                                 filter(!is_cortical) %>% 
                                 filter(cohort %in% c(cohorts[i],cohorts[j])), distribution=approximate(nresample=9999)))
    pw <- bind_rows(pw,tibble(param='log_knee',c1=cohorts[i],c2=cohorts[j],pvalue=as.numeric(pval)))
  }
}
pw_subcortical <- pw %>% mutate(is_cortical=0)

pw_adj <- bind_rows(pw_cortical,pw_subcortical) %>% 
  ungroup() %>%
  mutate(pvalue_adj = p.adjust(pvalue,method="fdr")) %>%
  write_tsv('data/per-cohort-pairwise-permutation-test-comparison')




