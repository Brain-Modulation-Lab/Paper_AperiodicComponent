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

UPDRS_laterality <- read_tsv('data/DBS_UPDRS_preop_laterality.tsv')

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
  left_join(UPDRS_laterality)

spearman_p <- function(x,y) cor.test(x,y,method = 'spearman')$p.value
spearman_r <- function(x,y) cor.test(x,y,method = 'spearman')$estimate


##### 

p_exponent_region1 <- DBS_aper_agg_region1 %>%
  group_by(region1) %>%
  summarise(p_off_axial = spearman_p(exponent,off_axial),
            p_off_left = spearman_p(exponent,off_left),
            p_off_right= spearman_p(exponent,off_right),
            p_on_axial = spearman_p(exponent,on_axial),
            p_on_left = spearman_p(exponent,on_left),
            p_on_right= spearman_p(exponent,on_right),
            r_off_axial = spearman_r(exponent,off_axial),
            r_off_left = spearman_r(exponent,off_left),
            r_off_right= spearman_r(exponent,off_right),
            r_on_axial = spearman_r(exponent,on_axial),
            r_on_left = spearman_r(exponent,on_left),
            r_on_right= spearman_r(exponent,on_right),
            n=n()) %>%
  filter(n>10)

p_exponent_region1_fdr <- p_exponent_region1 %>% 
  dplyr::select(!starts_with('r_')) %>%
  pivot_longer(!(region1|n),names_to="UPDRS_laterality",values_to="p_value") %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value,method="fdr")) %>%
  mutate(UPDRS_laterality=paste0(UPDRS_laterality,"_fdr")) %>%
  pivot_wider(id_cols=c(region1,n),names_from=UPDRS_laterality,values_from=p_value_fdr) %>%
  left_join(p_exponent_region1) 

p_exponent_region1_fdr %>%
  View()

# region1          n p_off…¹ p_off…² p_off…³ p_on_…⁴ p_on_…⁵ p_on_…⁶ p_off…⁷ p_off…⁸ p_off…⁹ p_on_…˟ p_on_…˟ p_on_…˟ r_off…˟ r_off…˟
# <chr>        <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
# 1 auditory-as…    21   0.745   0.848   0.699   0.699   0.745   0.745   0.390  0.803   0.324    0.335  0.544   0.479   0.198   0.0580
# 2 dorsolatera…    17   0.745   0.296   0.462   0.699   0.296   0.499   0.598  0.0424  0.0963   0.284  0.0311  0.114   0.138   0.497 
# 3 inferior-fr…    13   0.814   0.221   0.745   0.699   0.221   0.745   0.746  0.0103  0.434    0.318  0.0122  0.441   0.0996  0.681 
# 4 inferior-pa…    19   0.745   0.865   0.745   0.848   0.745   0.746   0.524  0.847   0.605    0.813  0.492   0.669   0.156  -0.0476
# 5 posterior-o…    20   0.965   0.746   0.745   0.745   0.699   0.745   0.965  0.662   0.559    0.600  0.292   0.553  -0.0106  0.104 
# 6 premotor        24   0.699   0.305   0.699   0.699   0.221   0.699   0.233  0.0508  0.247    0.254  0.0142  0.276   0.253   0.403 
# 7 sensorimotor    27   0.699   0.462   0.746   0.745   0.221   0.745   0.295  0.0937  0.634    0.419  0.0184  0.588   0.209   0.329 
# 8 temporo-par…    12   0.746   0.699   0.745   0.699   0.699   0.296   0.648  0.346   0.533    0.223  0.350   0.0432 -0.147  -0.299 
# # … with 4 more variables: r_off_right <dbl>, r_on_axial <dbl>, r_on_left <dbl>, r_on_right <dbl>, and abbreviated variable names
# #   ¹​p_off_axial_fdr, ²​p_off_left_fdr, ³​p_off_right_fdr, ⁴​p_on_axial_fdr, ⁵​p_on_left_fdr, ⁶​p_on_right_fdr, ⁷​p_off_axial,
# #   ⁸​p_off_left, ⁹​p_off_right, ˟​p_on_axial, ˟​p_on_left, ˟​p_on_right, ˟​r_off_axial, ˟​r_off_left
# # ℹ Use `colnames()` to see all variable names
#   
#knee 





#####

UPDRS_latsub <- read_tsv('data/DBS_UPDRS_preop_lat_sub.tsv')

DBS_aper_agg_region1 <- DBS_aper_coord_dx %>%
  left_join(MMP1_regions) %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  filter(dbs_target=='STN') %>%
  filter(rsquared > 0.98) %>%
  group_by(subject,region1) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent),n=n()) %>%
  left_join(UPDRS_latsub)


p_exponent_region1 <- DBS_aper_agg_region1 %>%
  group_by(region1) %>%
  summarise(p_off_axial = spearman_p(exponent,off_axial),
            p_off_left_bradykinesia = spearman_p(exponent,off_left_bradykinesia),
            p_off_left_rigity = spearman_p(exponent,off_left_rigidity),
            p_off_left_tremor = spearman_p(exponent,off_left_tremor),
            p_off_right_bradykinesia = spearman_p(exponent,off_right_bradykinesia),
            p_off_right_rigity = spearman_p(exponent,off_right_rigidity),
            p_off_right_tremor = spearman_p(exponent,off_right_tremor),
            
            p_on_axial = spearman_p(exponent,on_axial),
            p_on_left_bradykinesia = spearman_p(exponent,on_left_bradykinesia),
            p_on_left_rigity = spearman_p(exponent,on_left_rigidity),
            p_on_left_tremor = spearman_p(exponent,on_left_tremor),
            p_on_right_bradykinesia = spearman_p(exponent,on_right_bradykinesia),
            p_on_right_rigity = spearman_p(exponent,on_right_rigidity),
            p_on_right_tremor = spearman_p(exponent,on_right_tremor),
            
            r_off_axial = spearman_r(exponent,off_axial),
            r_off_left_bradykinesia = spearman_r(exponent,off_left_bradykinesia),
            r_off_left_rigity = spearman_r(exponent,off_left_rigidity),
            r_off_left_tremor = spearman_r(exponent,off_left_tremor),
            r_off_right_bradykinesia = spearman_r(exponent,off_right_bradykinesia),
            r_off_right_rigity = spearman_r(exponent,off_right_rigidity),
            r_off_right_tremor = spearman_r(exponent,off_right_tremor),
            
            r_on_axial = spearman_r(exponent,on_axial),
            r_on_left_bradykinesia = spearman_r(exponent,on_left_bradykinesia),
            r_on_left_rigity = spearman_r(exponent,on_left_rigidity),
            r_on_left_tremor = spearman_r(exponent,on_left_tremor),
            r_on_right_bradykinesia = spearman_r(exponent,on_right_bradykinesia),
            r_on_right_rigity = spearman_r(exponent,on_right_rigidity),
            r_on_right_tremor = spearman_r(exponent,on_right_tremor),
            
            n=n()) %>%
  filter(n>10)

p_exponent_region1_fdr <- p_exponent_region1 %>% 
  dplyr::select(!starts_with('r_')) %>%
  pivot_longer(!(region1|n),names_to="UPDRS_lat_sub",values_to="p_value") %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value,method="fdr")) %>%
  mutate(UPDRS_lat_sub=paste0(UPDRS_lat_sub,"_fdr")) %>%
  pivot_wider(id_cols=c(region1,n),names_from=UPDRS_lat_sub,values_from=p_value_fdr) %>%
  left_join(p_exponent_region1) 

p_exponent_region1_fdr %>%
  View()

# # A tibble: 8 × 44
# region1          n p_off…¹ p_off…² p_off…³ p_off…⁴ p_off…⁵ p_off…⁶ p_off…⁷ p_on_…⁸ p_on_…⁹ p_on_…˟ p_on_…˟ p_on_…˟ p_on_…˟ p_on_…˟
# <chr>        <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
# 1 auditory-as…    21   0.801  0.987    0.973   0.802   0.987   0.987   0.631   0.728  0.987    0.987  0.318    0.987   0.973   0.335
# 2 dorsolatera…    17   0.897  0.302    0.545   0.728   0.335   0.813   0.728   0.712  0.172    0.795  0.369    0.158   0.897   0.987
# 3 inferior-fr…    13   0.973  0.0282   0.397   0.973   0.431   0.987   0.973   0.728  0.0490   0.359  0.631    0.312   0.930   0.631
# 4 inferior-pa…    19   0.890  0.973    0.802   0.987   0.897   0.987   0.881   0.987  0.987    0.931  0.512    0.801   0.987   0.930
# 5 posterior-o…    20   0.987  0.987    0.987   0.634   0.890   0.987   0.369   0.897  0.987    0.987  0.172    0.987   0.965   0.397
# 6 premotor        24   0.669  0.369    0.897   0.335   0.712   0.987   0.369   0.682  0.369    0.819  0.0128   0.669   0.973   0.512
# 7 sensorimotor    27   0.712  0.397    0.868   0.506   0.802   0.931   0.499   0.802  0.379    0.728  0.158    0.682   0.897   0.890
# 8 temporo-par…    12   0.930  0.897    0.712   0.987   0.881   0.397   0.987   0.669  0.682    0.865  0.712    0.335   0.359   0.951
# # … with 28 more variables: p_off_axial <dbl>, p_off_left_bradykinesia <dbl>, p_off_left_rigity <dbl>, p_off_left_tremor <dbl>,
# #   p_off_right_bradykinesia <dbl>, p_off_right_rigity <dbl>, p_off_right_tremor <dbl>, p_on_axial <dbl>,
# #   p_on_left_bradykinesia <dbl>, p_on_left_rigity <dbl>, p_on_left_tremor <dbl>, p_on_right_bradykinesia <dbl>,
# #   p_on_right_rigity <dbl>, p_on_right_tremor <dbl>, r_off_axial <dbl>, r_off_left_bradykinesia <dbl>, r_off_left_rigity <dbl>,
# #   r_off_left_tremor <dbl>, r_off_right_bradykinesia <dbl>, r_off_right_rigity <dbl>, r_off_right_tremor <dbl>, r_on_axial <dbl>,
# #   r_on_left_bradykinesia <dbl>, r_on_left_rigity <dbl>, r_on_left_tremor <dbl>, r_on_right_bradykinesia <dbl>,
# #   r_on_right_rigity <dbl>, r_on_right_tremor <dbl>, and abbreviated variable names ¹​p_off_axial_fdr, …
# # ℹ Use `colnames()` to see all variable names

UPDRS_latsub %>%
  pivot_longer(cols=!subject) %>%
  separate(name,c('medstate','laterality','subscore')) %>%
  filter(!is.na(subscore)) %>%
  pivot_wider(id_cols=c(subject,medstate,subscore),names_from=laterality,values_from='value') %>%
  ggplot() +
  aes(x=right,y=left,color=medstate) +
  geom_point()+
  facet_wrap(~subscore) +
  geom_abline(slope=1,intercept=0)
ggsave('fig/B12c_01_UPDRS-score-lateralization.png')




