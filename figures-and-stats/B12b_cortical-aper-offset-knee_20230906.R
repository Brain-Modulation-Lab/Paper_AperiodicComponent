library(tidyverse)
library(coin)
library(multcomp)
library(lmerTest)
library(ggExtra)
library(ggradar)


theme_set(theme_bw())

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')

subjects %>%
  group_by(dx,dbs_target) %>%
  summarise(n=n())

artifacts <- read_tsv('data/artifacts.txt') %>%
  fill(subject,session_id,type) %>%
  mutate(subject=str_replace(subject,fixed('-'),'')) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_'))

# artifacts %>% write_tsv('data/artifacts.txt')

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
  mutate(electrode=tolower(electrode)) %>%
  left_join(subjects %>% dplyr::select(subject,dx,dbs_target)) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel)

DBS_aper_coord_dx <- read_tsv('data/A01_DBS_aper_coord_dx.tsv') %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel)

MMP1_regions <- read_tsv('data/MMP1-regions2.tsv')
fs_anatomy_lobules <- read_tsv('data/fs_anatomy_lobules.tsv')
fs_anatomy_regions<- read_tsv('data/fs_anatomy_regions.tsv')

UPDRS_subscores <- read_tsv('data/DBS_UPDRS_preop_subscores.tsv')

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG


##############
# 
# DBS_aper_coord_dx %>%
#   filter(electrode_type=='ecog') %>%
#   #filter(dx=='PD') %>%
#   with(tibble(name=setdiff(unique(HCPMMP1_label_1),MMP1_regions$HCPMMP1_label_1))) %>%
#   View()

DBS_aper_agg_region1 <- DBS_aper_coord_dx %>%
  left_join(MMP1_regions) %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  filter(dbs_target=='STN') %>%
  filter(rsquared > 0.98) %>%
  group_by(subject,region1) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent),n=n()) %>%
  left_join(UPDRS_subscores)

spearman_p <- function(x,y) cor.test(x,y,method = 'spearman')$p.value
spearman_r <- function(x,y) cor.test(x,y,method = 'spearman')$estimate


##### Offset

p_offset_region1 <- DBS_aper_agg_region1 %>%
  group_by(region1) %>%
  summarise(p_off_total = spearman_p(offset,off_total),
            p_on_total = spearman_p(offset,on_total),
            r_off_total = spearman_r(offset,off_total),
            r_on_total = spearman_r(offset,on_total),
            n=n()) %>%
  filter(n>10)

p_offset_region1_fdr <- p_offset_region1 %>% 
  dplyr::select(!starts_with('r_')) %>%
  pivot_longer(!(region1|n),names_to="UPDRS_score",values_to="p_value") %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value,method="fdr")) %>%
  pivot_wider(id_cols=c(region1,n),names_from=UPDRS_score,values_from=p_value_fdr) %>%
  rename(p_off_total_fdr=p_off_total,p_on_total_fdr=p_on_total) %>%
  left_join(p_offset_region1) 

p_offset_region1_fdr %>%
  View()

# # A tibble: 8 × 8
# region1                        n p_off_total_fdr p_on_total_fdr p_off_total p_on_total r_off_total r_on_total
# <chr>                      <int>           <dbl>          <dbl>       <dbl>      <dbl>       <dbl>      <dbl>
# 1 auditory-associative          21           0.559          0.559      0.323      0.308       0.227      0.234 
# 2 dorsolateral-prefrontal       16           0.559          0.185      0.350      0.0232      0.250      0.563 
# 3 inferior-frontal              12           0.714          0.343      0.704      0.0945      0.123      0.504 
# 4 inferior-parietal             18           0.714          0.714      0.671      0.586      -0.107      0.138 
# 5 posterior-opercular           20           0.714          0.714      0.714      0.684       0.0873     0.0972
# 6 premotor                      23           0.195          0.185      0.0365     0.0196      0.438      0.483 
# 7 sensorimotor                  26           0.549          0.343      0.206      0.107       0.257      0.323 
# 8 temporo-parietal-occipital    12           0.714          0.559      0.619      0.269       0.161      0.347 

  
#knee 

p_log_knee_region1 <- DBS_aper_agg_region1 %>%
  group_by(region1) %>%
  summarise(p_off_total = spearman_p(log_knee,off_total),
            p_on_total = spearman_p(log_knee,on_total),
            r_off_total = spearman_r(log_knee,off_total),
            r_on_total = spearman_r(log_knee,on_total),
            n=n()) %>%
  filter(n>10)

p_log_knee_region1_fdr <- p_log_knee_region1 %>% 
  dplyr::select(!starts_with('r_')) %>%
  pivot_longer(!(region1|n),names_to="UPDRS_score",values_to="p_value") %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value,method="fdr")) %>%
  pivot_wider(id_cols=c(region1,n),names_from=UPDRS_score,values_from=p_value_fdr) %>%
  rename(p_off_total_fdr=p_off_total,p_on_total_fdr=p_on_total) %>%
  left_join(p_log_knee_region1) 

p_log_knee_region1_fdr %>%
  View()

# # A tibble: 8 × 8
# region1                        n p_off_total_fdr p_on_total_fdr p_off_total p_on_total r_off_total r_on_total
# <chr>                      <int>           <dbl>          <dbl>       <dbl>      <dbl>       <dbl>      <dbl>
# 1 auditory-associative          21           0.882          0.882      0.735      0.827      -0.0786    -0.0507
# 2 dorsolateral-prefrontal       16           0.882          0.882      0.824      0.590      -0.0604     0.146 
# 3 inferior-frontal              12           0.449          0.449      0.0842     0.0682      0.518      0.543 
# 4 inferior-parietal             18           0.878          0.882      0.440      0.671      -0.194     -0.108 
# 5 posterior-opercular           20           0.878          0.692      0.383      0.173       0.206      0.317 
# 6 premotor                      23           0.957          0.878      0.957      0.469       0.0119     0.159 
# 7 sensorimotor                  26           0.882          0.878      0.752      0.494       0.0650     0.140 
# 8 temporo-parietal-occipital    12           0.878          0.449      0.457      0.0794     -0.238     -0.525 

