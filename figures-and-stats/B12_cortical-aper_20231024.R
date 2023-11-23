library(tidyverse)
library(coin)
library(multcomp)
library(lmerTest)
library(ggExtra)
library(ggradar)


theme_set(theme_bw())

#setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')
setwd('~/Dropbox/Lab-BML/Expan/2021-11-16-FOOOF-figures')



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

p_exp_region1 <- DBS_aper_agg_region1 %>%
  group_by(region1) %>%
  summarise(p_off_total = spearman_p(exponent,off_total),
            p_on_total = spearman_p(exponent,on_total),
            r_off_total = spearman_r(exponent,off_total),
            r_on_total = spearman_r(exponent,on_total),
            n=n()) %>%
  filter(n>10)

p_exp_region1_fdr <- p_exp_region1 %>% 
  dplyr::select(!starts_with('r_')) %>%
  pivot_longer(!(region1|n),names_to="UPDRS_score",values_to="p_value") %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value,method="fdr")) %>%
  pivot_wider(id_cols=c(region1,n),names_from=UPDRS_score,values_from=p_value_fdr) %>%
  rename(p_off_total_fdr=p_off_total,p_on_total_fdr=p_on_total) %>%
  left_join(p_exp_region1) 

p_exp_region1_fdr %>%
  View()
#   region1                        n p_off_total_fdr p_on_total_fdr p_off_total p_on_total r_off_total r_on_total
# 1 auditory-associative          21          0.789          0.489       0.691     0.367        0.0923     0.208 
# 2 dorsolateral-prefrontal       16          0.429          0.0174      0.241     0.00217      0.311      0.707 
# 3 inferior-frontal              12          0.0790         0.0161      0.0247    0.00100      0.641      0.823 
# 4 inferior-parietal             18          0.547          0.961       0.445     0.961       -0.192      0.0124
# 5 posterior-opercular           20          0.941          0.489       0.882     0.362        0.0354     0.215 
# 6 premotor                      23          0.226          0.0507      0.0864    0.00951      0.365      0.529 
# 7 sensorimotor                  26          0.226          0.0690      0.113     0.0173       0.318      0.463 
# 8 temporo-parietal-occipital    12          0.489          0.226       0.331     0.0999      -0.308     -0.497 

p_exp_subscores <- DBS_aper_agg_region1 %>%
  filter(region1 %in% c("dorsolateral-prefrontal","inferior-frontal")) %>%
  group_by(region1) %>%
  summarise(
    p_on_axial = spearman_p(exponent,on_axial),
    p_on_bradykinesia = spearman_p(exponent,on_bradykinesia),            
    p_on_rigidity = spearman_p(exponent,on_rigidity),     
    p_on_tremor = spearman_p(exponent,on_tremor),
    r_on_axial = spearman_r(exponent,on_axial),
    r_on_bradykinesia = spearman_r(exponent,on_bradykinesia),            
    r_on_rigidity = spearman_r(exponent,on_rigidity),     
    r_on_tremor = spearman_r(exponent,on_tremor),
    
    p_off_axial = spearman_p(exponent,off_axial),
    p_off_bradykinesia = spearman_p(exponent,off_bradykinesia),            
    p_off_rigidity = spearman_p(exponent,off_rigidity),     
    p_off_tremor = spearman_p(exponent,off_tremor),
    r_off_axial = spearman_r(exponent,off_axial),
    r_off_bradykinesia = spearman_r(exponent,off_bradykinesia),            
    r_off_rigidity = spearman_r(exponent,off_rigidity),     
    r_off_tremor = spearman_r(exponent,off_tremor),
    
    n=n()
    ) %>%
  filter(n>10)

p_exp_subscores_fdr <- p_exp_subscores %>% 
  dplyr::select(!starts_with('r_')) %>%
  pivot_longer(!(region1|n),names_to="UPDRS_subscore",values_to="p_value") %>%
  ungroup() %>%
  mutate(p_value_fdr = p.adjust(p_value,method="fdr")) %>%
  mutate(UPDRS_subscore=paste0(UPDRS_subscore,"_fdr")) %>%
  pivot_wider(id_cols=c(region1,n),names_from=UPDRS_subscore,values_from=p_value_fdr) %>%
  left_join(p_exp_subscores) 

p_exp_subscores_fdr %>%
  View()

#   region1        n p_on_…¹ p_on_…² p_on_…³ p_on_…⁴ p_off…⁵ p_off…⁶ p_off…⁷ p_off…⁸ p_on_…⁹ p_on_…˟ p_on_…˟ p_on_…˟
#     <chr>      <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
# 1 dorsolate…    16   0.277 0.0327   0.277    0.277   0.593  0.166    0.277   0.277   0.184 4.55e-3  0.144    0.208
# 2 inferior-…    12   0.277 0.00796  0.0444   0.909   0.759  0.0327   0.277   0.931   0.150 4.98e-4  0.0111   0.852


#============= plotting

myTheme <- theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                              strip.background = element_blank())

p_exp_region1_fdr_label <- p_exp_region1_fdr %>%
  mutate(label_on=paste0('p=',signif(p_on_total_fdr,2),' ','rho=',signif(r_on_total,2))) %>%
  mutate(label_off=paste0('p=',signif(p_off_total_fdr,2),' ','rho=',signif(r_off_total,2)))
  

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_region1_fdr$region1) %>%
  ggplot() +
  aes(x=on_total,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,35)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_on,x=NULL,Y=NULL),
             data=p_exp_region1_fdr_label,
             x=35,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_11a_exponent-vs-UPDRS-on-per-region.pdf',width=7,height=4.5)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_region1_fdr$region1) %>%
  ggplot() +
  aes(x=off_total,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_off,x=NULL,Y=NULL),
             data=p_exp_region1_fdr_label,
             x=50,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_11b_exponent-vs-UPDRS-off-per-region.pdf',width=7,height=4.5)


# subscores
p_exp_subscores_fdr_labels <- p_exp_subscores_fdr %>%
  mutate(label_on_axial=paste0('p=',signif(p_on_axial_fdr,2),' ','rho=',signif(r_on_axial,2))) %>%
  mutate(label_on_bradykinesia=paste0('p=',signif(p_on_bradykinesia_fdr,2),' ','rho=',signif(r_on_bradykinesia,2))) %>%
  mutate(label_on_rigidity=paste0('p=',signif(p_on_rigidity_fdr,2),' ','rho=',signif(r_on_rigidity,2))) %>%
  mutate(label_on_tremor=paste0('p=',signif(p_on_tremor_fdr,2),' ','rho=',signif(r_on_tremor,2))) %>%
  mutate(label_off_axial=paste0('p=',signif(p_off_axial_fdr,2),' ','rho=',signif(r_off_axial,2))) %>%
  mutate(label_off_bradykinesia=paste0('p=',signif(p_off_bradykinesia_fdr,2),' ','rho=',signif(r_off_bradykinesia,2))) %>%
  mutate(label_off_rigidity=paste0('p=',signif(p_off_rigidity_fdr,2),' ','rho=',signif(r_off_rigidity,2))) %>%
  mutate(label_off_tremor=paste0('p=',signif(p_off_tremor_fdr,2),' ','rho=',signif(r_off_tremor,2)))

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=on_axial,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_on_axial,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=10,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-on-axial-per-region.pdf',width=4,height=3)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=on_bradykinesia,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,20)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_on_bradykinesia,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=20,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-on-bradykinesia-per-region.pdf',width=4,height=3)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=on_rigidity,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,12)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_on_rigidity,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=12,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-on-rigidity-per-region.pdf',width=4,height=3)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=on_tremor,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_on_tremor,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=10,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-on-tremor-per-region.pdf',width=4,height=3)



DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=off_axial,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,15)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_off_axial,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=15,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-off-axial-per-region.pdf',width=4,height=3)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=off_bradykinesia,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_off_bradykinesia,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=30,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-off-bradykinesia-per-region.pdf',width=4,height=3)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=off_rigidity,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,15)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_off_rigidity,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=15,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-off-rigidity-per-region.pdf',width=4,height=3)

DBS_aper_agg_region1 %>%
  filter(region1 %in% p_exp_subscores_fdr$region1) %>%
  ggplot() +
  aes(x=off_tremor,y=exponent)+
  geom_smooth(method='lm',size=0.5,color='black') +
  geom_point() +
  facet_wrap(~region1,scales='free') +
  scale_x_continuous(limits=c(0,15)) + scale_y_continuous(limits=c(2,4.5)) +
  geom_label(mapping=aes(label=label_off_tremor,x=NULL,Y=NULL),
             data=p_exp_subscores_fdr_labels,
             x=15,y=2,hjust=1,vjust=0,label.size = NA,size=3) +
  myTheme
ggsave('fig/B12_12a_exponent-vs-UPDRS-off-tremor-per-region.pdf',width=4,height=3)



radar_db <- p_exp_subscores_fdr %>%
  dplyr::select(region1,starts_with('r_')) %>%
  pivot_longer(starts_with('r_')) %>%
  separate(name,c(NA,'state','subscore')) %>%
  pivot_wider(id_cols=c(region1,state),names_from=subscore) 


for (ur in unique(radar_db$region1)){
  p1 <- radar_db %>%
  filter(region1==ur) %>%
  dplyr::select(-region1) %>%
  ggradar(fill=TRUE,fill.alpha=0.25,values.radar = c("0", ".5", "1"),
          group.line.width = 1,  group.point.size = 5) +
  theme_void()
  
  ggsave(paste0('fig/B12_13a_radar-plot_',ur,'.pdf'),plot=p1,width=4,height=4)
}


###### Table with number of "good" cortical electrodes per region per subject

DBS_aper_agg_region1_n <- DBS_aper_coord_dx %>%
  left_join(MMP1_regions) %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(rsquared > 0.98) %>%
  group_by(subject,session_id,dx,dbs_target,region1) %>%
  summarise(n=n()) %>%
  group_by(subject,dx,dbs_target,region1) %>%
  summarise(n=max(n)) %>%
  ungroup() %>%
  pivot_wider(id_cols=c(subject,dx,dbs_target),names_from=region1,values_from=n,values_fill=0) %>%
  left_join(UPDRS_subscores) %>%
  write_tsv("data/B12_summary-of-subjects.tsv")

setdiff(subjects$subject,DBS_aper_agg_region1_n$subject)
setdiff(subjects$subject,DBS_aper_coord_dx$subject)

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv')
setdiff(subjects$subject,DBS_aper$subject)



#### plotting parameters vs Dx

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG

DBS_aper_agg2 <- DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(!is.na(dbs_target)) %>%
  mutate(dx_dbs_target=str_c(dx,'_',dbs_target)) %>%
  group_by(dx_dbs_target,subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  mutate(dx_dbs_target=factor(dx_dbs_target,levels=c('PD_STN','PD_GPi','ET_VIM'),ordered=TRUE)) %>%
  write_tsv('fig/B12_22_dx_DBS_target.tsv')

setdiff(DBS_aper$subject,DBS_aper_agg2$subject)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,x=dx_dbs_target,fill=dx_dbs_target)+
  geom_dotplot(binaxis = "y", stackdir = "center",binwidth=0.05,binposition='all') +
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_y_continuous(limits=c(2.5,4.1))+
  scale_fill_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B12_22_exponent-vs-dx.pdf',width=1.8,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,x=dx_dbs_target,color=dx_dbs_target,fill=dx_dbs_target)+
  geom_boxplot(alpha=0,width=0.6) +
  geom_point(position=position_jitter(width=0.25),colour="black",shape=21)+
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_y_continuous(limits=c(2.5,4.1))+
  scale_color_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))+
  scale_fill_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B12_22b_exponent-vs-dx.pdf',width=1.8,height=3)





DBS_aper_agg2 %>%
  group_by(dx_dbs_target) %>%
  summarise(n=n())

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=log_knee,x=dx_dbs_target,color=dx_dbs_target)+
  geom_boxplot(alpha=0,width=0.6) +
  geom_jitter(width=0.25)+
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_color_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B12_23b_log-knee-vs-dx.pdf',width=1.8,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=offset,x=dx_dbs_target,color=dx_dbs_target)+
  geom_boxplot(alpha=0,width=0.6) +
  geom_jitter(width=0.25)+
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_color_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B12_23b_offset-vs-dx.pdf',width=1.8,height=3)

#---------

summary(lm(exponent ~ dx_dbs_target, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      3.16327    0.06565  48.183   <2e-16 ***
# dx_dbs_target.L -0.12107    0.07280  -1.663   0.1033    
# dx_dbs_target.Q  0.28332    0.14339   1.976   0.0543 .  

oneway_test(exponent ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2, 
            distribution=approximate(nresample=999))
# Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by
# factor(dx_dbs_target, ordered = FALSE) (PD_STN, PD_GPi, ET_VIM)
# chi-squared = 6.3582, p-value = 0.03103

oneway_test(exponent ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2 %>% filter(dx_dbs_target %in% c("PD_STN","ET_VIM")), 
            distribution=approximate(nresample=999))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by factor(dx_dbs_target, ordered = FALSE) (PD_STN, ET_VIM)
# Z = 1.6246, p-value = 0.0961
# alternative hypothesis: true mu is not equal to 0

summary(lm(log_knee ~ dx_dbs_target, data=DBS_aper_agg2))
# Call:
#   lm(formula = log_knee ~ dx_dbs_target, data = DBS_aper_agg2)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.36380 -0.03220  0.01856  0.07956  0.16545 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.155907   0.025182  45.903  < 2e-16 ***
#   dx_dbs_target.L -0.005973   0.027926  -0.214  0.83161    
# dx_dbs_target.Q  0.167852   0.054999   3.052  0.00381 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1288 on 45 degrees of freedom
# Multiple R-squared:  0.1737,	Adjusted R-squared:  0.1369 
# F-statistic: 4.728 on 2 and 45 DF,  p-value: 0.01368

oneway_test(log_knee ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2, 
            distribution=approximate(nresample=999))
# Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by
# factor(dx_dbs_target, ordered = FALSE) (PD_STN, PD_GPi, ET_VIM)
# chi-squared = 8.1616, p-value = 0.01802

oneway_test(log_knee ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2 %>% filter(dx_dbs_target %in% c("PD_STN","ET_VIM")), 
            distribution=approximate(nresample=999))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by factor(dx_dbs_target, ordered = FALSE) (PD_STN, ET_VIM)
# Z = 0.22837, p-value = 0.8378
# alternative hypothesis: true mu is not equal to 0

summary(lm(offset ~ dx_dbs_target, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       1.7085     0.1070  15.963   <2e-16 ***
#   dx_dbs_target.L  -0.1192     0.1187  -1.004    0.321    
# dx_dbs_target.Q  -0.1374     0.2338  -0.588    0.560    

oneway_test(offset ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2, 
            distribution=approximate(nresample=999))
# Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  offset by
# factor(dx_dbs_target, ordered = FALSE) (PD_STN, PD_GPi, ET_VIM)
# chi-squared = 1.3146, p-value = 0.5185

#=======================



###### Correlation with age


DBS_aper_agg24 <- DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(!is.na(dbs_target)) %>%
  mutate(dx_dbs_target=str_c(dx,'_',dbs_target)) %>%
  dplyr::filter(dx_dbs_target %in% c('PD_STN', 'ET_VIM')) %>%
  group_by(dx_dbs_target,subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects %>% dplyr::select(subject,dbs_age)) %>%
  write_tsv('fig/B12_24_aper-vs-age.tsv')

DBS_aper_agg24 %>%
  ggplot() +
  aes(x=dbs_age,y=exponent) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('fig/B12_24_exponent-vs-age.pdf',width=5,height=4)

DBS_aper_agg24 %>%
  with(cor.test(exponent,dbs_age,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  exponent and dbs_age
# S = 11583, p-value = 0.2326
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1837046 

DBS_aper_agg24 %>%
  ggplot() +
  aes(x=dbs_age,y=log_knee) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('fig/B12_24_log-knee-vs-age.pdf',width=5,height=4)

DBS_aper_agg24 %>%
  with(cor.test(log_knee,dbs_age,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  log_knee and dbs_age
# S = 14985, p-value = 0.7179
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.05603591 

DBS_aper_agg24 %>%
  ggplot() +
  aes(x=dbs_age,y=offset) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('fig/B12_24_offset-vs-age.pdf',width=5,height=4)

DBS_aper_agg24 %>%
  with(cor.test(offset,dbs_age,method = 'spearman'))

# Spearman's rank correlation rho
# 
# data:  offset and dbs_age
# S = 12121, p-value = 0.345
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1458063 


########### doing correlation vs anatomy

coverage <- DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(!is.na(dbs_target)) %>%
  mutate(dx_dbs_target=str_c(dx,'_',dbs_target)) %>%
  filter(dx_dbs_target%in%c('PD_STN','ET_VIM')) %>%
  group_by(subject, HCPMMP1_label_1) %>%
  tally(name='n_e') %>%
  group_by(HCPMMP1_label_1) %>%
  summarise(n_e=sum(n_e),n_s=n()) %>%
  arrange(desc(n_s),desc(n_e))

roi <- coverage %>%
  filter(n_s >= 10)

aper_ecog_roi <-  DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(!is.na(dbs_target)) %>%
  mutate(dx_dbs_target=str_c(dx,'_',dbs_target)) %>%
  filter(dx_dbs_target%in%c('PD_STN','ET_VIM')) %>%
  filter(HCPMMP1_label_1 %in% roi$HCPMMP1_label_1)

#checking that all subjects have sufficient regions
aper_ecog_roi %>%
  group_by(subject, HCPMMP1_label_1) %>%
  tally(name='n_e') %>%
  group_by(subject) %>%
  summarise(n_e=sum(n_e),n_r=n()) %>%
  arrange(desc(n_r),desc(n_e)) %>%
  View()

aper_ecog_roi_agg <- aper_ecog_roi %>%
  filter(log_knee > 1) %>%
  group_by(subject,HCPMMP1_label_1) %>%
  summarise(exponent = median(exponent), log_knee = median(log_knee), offset = median(offset),
            x=mean(mni_nonlinear_x),y=mean(mni_nonlinear_y),z=mean(mni_nonlinear_z)) %>%
  group_by(subject) %>% 
  mutate(exponent_s = exponent - median(exponent), log_knee_s = log_knee - median(log_knee), offset_s = offset - median(offset)) %>%
  ungroup() %>%
  mutate(region = factor(HCPMMP1_label_1,levels=MMP1_regions$HCPMMP1_label_1)) %>%
  left_join(MMP1_regions) %>%
  mutate(region1 = factor(region1,levels=unique(MMP1_regions$region1))) %>%
  arrange(region)

grand_median_for_plotting <- aper_ecog_roi %>%
  filter(log_knee > 1) %>%
  group_by(subject,HCPMMP1_label_1) %>%
  summarise(exponent = median(exponent), log_knee = median(log_knee), offset = median(offset),
            x=mean(mni_nonlinear_x),y=mean(mni_nonlinear_y),z=mean(mni_nonlinear_z)) %>%
  group_by(subject) %>% 
  summarise(median_exponent_s = median(exponent), median_log_knee_s = median(log_knee), median_offset_s = median(offset)) %>%
  ungroup() %>%
  summarise(median_exponent_s = median(median_exponent_s), median_log_knee_s = median(median_log_knee_s), median_offset_s = median(median_offset_s)) 
# median_exponent_s median_log_knee_s median_offset_s
# <dbl>             <dbl>           <dbl>
#   1              3.36              1.25            1.72


aper_ecog_roi_agg %>%
  ggplot() + 
  aes(x=region,y=exponent,color=subject,group=subject)+
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))

aper_ecog_roi_agg %>%
  ggplot() + 
  aes(x=region,y=exponent_s,color=region1)+
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))

aper_ecog_roi_agg %>%
  filter(log_knee > 1.1) %>%
  ggplot() + 
  aes(x=region,y=log_knee,color=subject,group=subject)+
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))

aper_ecog_roi_agg %>%
  ggplot() + 
  aes(x=region,y=log_knee_s,color=region1)+
  geom_boxplot() + 
  scale_color_brewer(type='qual',palette='Set2') +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  theme(legend.position='top',legend.title = element_blank(),
        axis.title.x=element_blank(),legend.text=element_text(size=6)) +
  labs(y="Log10(Knee frequency)\n region effect")
ggsave('fig/B12_28_aper_log-knee_vs_anatomy.pdf',width=5,height=4)

aper_ecog_roi_agg %>%
  group_by(region,region1) %>%
  summarise(mean_log_knee_s = mean(log_knee_s,na.rm=TRUE), se_log_knee_s=sd(log_knee_s,na.rm=TRUE)/sqrt(n())) %>%
  ggplot() + 
  aes(x=region,y=mean_log_knee_s,color=region1,ymax=mean_log_knee_s+se_log_knee_s,ymin=mean_log_knee_s-se_log_knee_s)+
  geom_errorbar() +
  geom_point()+
  scale_color_brewer(type='qual',palette='Set2') +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  theme(legend.position='top',legend.title = element_blank(),
        axis.title.x=element_blank(),legend.text=element_text(size=6)) +
  labs(y="Log10(Knee frequency)\n region effect")
ggsave('fig/B12_28b_aper_log-knee_vs_anatomy.pdf',width=5,height=4)

param_stat_by_MMP1 <- aper_ecog_roi_agg %>%
  group_by(region,region1) %>%
  summarise(n_s=length(unique(subject)), 
            mean_log_knee_s = mean(log_knee_s,na.rm=TRUE), se_log_knee_s=sd(log_knee_s,na.rm=TRUE)/sqrt(n()),
            mean_exponent_s = mean(exponent_s,na.rm=TRUE), se_exponent_s=sd(exponent_s,na.rm=TRUE)/sqrt(n()),
            mean_offset_s = mean(offset_s,na.rm=TRUE), se_offset_s=sd(offset_s,na.rm=TRUE)/sqrt(n())) #%>%
param_stat_by_MMP1 <- aper_ecog_roi_agg %>%
  group_by(region,region1,subject) %>%
  summarise(median_log_knee = median(log_knee,na.rm=TRUE),median_exponent = median(exponent,na.rm=TRUE),median_offset = median(offset,na.rm=TRUE)) %>%
  group_by(region,region1) %>%
  summarise(    
    mean_log_knee = mean(median_log_knee,na.rm=TRUE), se_log_knee=sd(median_log_knee,na.rm=TRUE)/sqrt(n()),
    mean_exponent = mean(median_exponent,na.rm=TRUE), se_exponent=sd(median_exponent,na.rm=TRUE)/sqrt(n()),
    mean_offset = mean(median_offset,na.rm=TRUE), se_offset=sd(median_offset,na.rm=TRUE)/sqrt(n())) %>%
  mutate(ci_up_knee_frequency = 10^(mean_log_knee+se_log_knee), ci_low_knee_frequency = 10^(mean_log_knee-se_log_knee),
         mean_knee_frequency = 10^mean_log_knee, se_knee_frequency=(ci_up_knee_frequency-ci_low_knee_frequency)/2) %>%
  left_join(param_stat_by_MMP1) %>%
  dplyr::select(region1,region,n_s,everything()) %>%
  write_tsv('data/B12_aper_param_summary_by_MMP1.tsv')

me_mod <- lmer(log_knee ~ region + (1 | subject), 
               data=aper_ecog_roi_agg %>% within(region<-relevel(region,'6v')))
summary(me_mod)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log_knee ~ region + (1 | subject)
#    Data: aper_ecog_roi_agg %>% within(region <- relevel(region, "6v"))
# 
# REML criterion at convergence: -807.7
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.8198 -0.5802  0.0489  0.5377  3.8646 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  subject  (Intercept) 0.005468 0.07395 
#  Residual             0.005389 0.07341 
# Number of obs: 425, groups:  subject, 44
# 
# Fixed effects:
#                Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)   1.257e+00  1.627e-02  1.373e+02  77.235   <2e-16 ***
# regionp9-46v -5.795e-03  2.682e-02  3.653e+02  -0.216   0.8290    
# region8C     -1.613e-03  1.966e-02  3.637e+02  -0.082   0.9347    
# region8Av    -2.379e-02  2.057e-02  3.637e+02  -1.156   0.2483    
# regionIFSp   -3.391e-03  2.496e-02  3.652e+02  -0.136   0.8920    
# region44     -8.259e-03  2.231e-02  3.644e+02  -0.370   0.7115    
# region6r     -3.876e-02  2.181e-02  3.640e+02  -1.777   0.0764 .  
# region55b     8.173e-03  1.912e-02  3.635e+02   0.428   0.6692    
# region4      -1.616e-03  2.177e-02  3.637e+02  -0.074   0.9409    
# region3a     -1.644e-02  2.341e-02  3.641e+02  -0.702   0.4829    
# region3b      2.978e-02  2.002e-02  3.653e+02   1.488   0.1377    
# region1       3.891e-02  1.654e-02  3.620e+02   2.353   0.0192 *  
# region43      9.708e-04  2.092e-02  3.635e+02   0.046   0.9630    
# regionOP4    -2.786e-02  1.782e-02  3.620e+02  -1.564   0.1188    
# regionPFop   -1.925e-02  1.892e-02  3.640e+02  -1.017   0.3097    
# regionPF     -5.794e-03  2.206e-02  3.680e+02  -0.263   0.7930    
# regionPSL    -4.406e-02  2.822e-02  3.676e+02  -1.561   0.1194    
# regionSTV    -9.590e-03  2.602e-02  3.678e+02  -0.368   0.7127    
# regionA5     -4.078e-02  1.996e-02  3.640e+02  -2.044   0.0417 *  
# regionA4     -6.524e-03  1.732e-02  3.613e+02  -0.377   0.7067    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


amod <- aov(log_knee_s ~ region, data=aper_ecog_roi_agg)
amod
summary(amod)
tuk <- glht(amod, linfct = mcp(region = "Tukey"))
summary(tuk) 

# 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: aov(formula = log_knee_s ~ region, data = aper_ecog_roi_agg)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)   
# 8C - p9-46v == 0    1.325e-03  2.745e-02   0.048   1.0000   
# 8Av - p9-46v == 0  -2.053e-02  2.807e-02  -0.731   1.0000   
# IFSp - p9-46v == 0 -2.421e-04  3.103e-02  -0.008   1.0000   
# 44 - p9-46v == 0    1.053e-03  2.921e-02   0.036   1.0000   
# 6r - p9-46v == 0   -3.222e-02  2.888e-02  -1.116   0.9999   
# 6v - p9-46v == 0    2.793e-03  2.569e-02   0.109   1.0000   
# 55b - p9-46v == 0   1.041e-02  2.711e-02   0.384   1.0000   
# 4 - p9-46v == 0    -4.000e-03  2.888e-02  -0.138   1.0000   
# 3a - p9-46v == 0   -1.675e-02  3.000e-02  -0.558   1.0000   
# 3b - p9-46v == 0    2.801e-02  2.764e-02   1.014   1.0000   
# 1 - p9-46v == 0     4.408e-02  2.556e-02   1.725   0.9723   
# 43 - p9-46v == 0   -2.302e-03  2.831e-02  -0.081   1.0000   
# OP4 - p9-46v == 0  -2.514e-02  2.635e-02  -0.954   1.0000   
# PFop - p9-46v == 0 -1.455e-02  2.697e-02  -0.540   1.0000   
# PF - p9-46v == 0   -4.853e-03  2.888e-02  -0.168   1.0000   
# PSL - p9-46v == 0  -3.780e-02  3.330e-02  -1.135   0.9999   
# STV - p9-46v == 0  -1.439e-02  3.166e-02  -0.454   1.0000   
# A5 - p9-46v == 0   -3.961e-02  2.764e-02  -1.433   0.9966   
# A4 - p9-46v == 0   -3.973e-03  2.607e-02  -0.152   1.0000   
# 8Av - 8C == 0      -2.185e-02  2.216e-02  -0.986   1.0000   
# IFSp - 8C == 0     -1.567e-03  2.581e-02  -0.061   1.0000   
# 44 - 8C == 0       -2.720e-04  2.359e-02  -0.012   1.0000   
# 6r - 8C == 0       -3.354e-02  2.318e-02  -1.447   0.9962   
# 6v - 8C == 0        1.468e-03  1.905e-02   0.077   1.0000   
# 55b - 8C == 0       9.085e-03  2.094e-02   0.434   1.0000   
# 4 - 8C == 0        -5.325e-03  2.318e-02  -0.230   1.0000   
# 3a - 8C == 0       -1.808e-02  2.457e-02  -0.736   1.0000   
# 3b - 8C == 0        2.669e-02  2.161e-02   1.235   0.9995   
# 1 - 8C == 0         4.275e-02  1.888e-02   2.265   0.7437   
# 43 - 8C == 0       -3.627e-03  2.247e-02  -0.161   1.0000   
# OP4 - 8C == 0      -2.647e-02  1.994e-02  -1.327   0.9988   
# PFop - 8C == 0     -1.587e-02  2.074e-02  -0.765   1.0000   
# PF - 8C == 0       -6.178e-03  2.318e-02  -0.267   1.0000   
# PSL - 8C == 0      -3.912e-02  2.849e-02  -1.373   0.9981   
# STV - 8C == 0      -1.571e-02  2.657e-02  -0.591   1.0000   
# A5 - 8C == 0       -4.094e-02  2.161e-02  -1.894   0.9323   
# A4 - 8C == 0       -5.298e-03  1.956e-02  -0.271   1.0000   
# IFSp - 8Av == 0     2.028e-02  2.646e-02   0.767   1.0000   
# 44 - 8Av == 0       2.158e-02  2.431e-02   0.888   1.0000   
# 6r - 8Av == 0      -1.169e-02  2.391e-02  -0.489   1.0000   
# 6v - 8Av == 0       2.332e-02  1.993e-02   1.170   0.9998   
# 55b - 8Av == 0      3.094e-02  2.174e-02   1.423   0.9969   
# 4 - 8Av == 0        1.653e-02  2.391e-02   0.691   1.0000   
# 3a - 8Av == 0       3.774e-03  2.525e-02   0.149   1.0000   
# 3b - 8Av == 0       4.854e-02  2.239e-02   2.168   0.8074   
# 1 - 8Av == 0        6.460e-02  1.977e-02   3.269   0.1148   
# 43 - 8Av == 0       1.822e-02  2.322e-02   0.785   1.0000   
# OP4 - 8Av == 0     -4.616e-03  2.078e-02  -0.222   1.0000   
# PFop - 8Av == 0     5.977e-03  2.155e-02   0.277   1.0000   
# PF - 8Av == 0       1.567e-02  2.391e-02   0.656   1.0000   
# PSL - 8Av == 0     -1.727e-02  2.909e-02  -0.594   1.0000   
# STV - 8Av == 0      6.137e-03  2.720e-02   0.226   1.0000   
# A5 - 8Av == 0      -1.909e-02  2.239e-02  -0.853   1.0000   
# A4 - 8Av == 0       1.655e-02  2.042e-02   0.811   1.0000   
# 44 - IFSp == 0      1.295e-03  2.767e-02   0.047   1.0000   
# 6r - IFSp == 0     -3.198e-02  2.732e-02  -1.170   0.9998   
# 6v - IFSp == 0      3.035e-03  2.392e-02   0.127   1.0000   
# 55b - IFSp == 0     1.065e-02  2.545e-02   0.419   1.0000   
# 4 - IFSp == 0      -3.758e-03  2.732e-02  -0.138   1.0000   
# 3a - IFSp == 0     -1.651e-02  2.851e-02  -0.579   1.0000   
# 3b - IFSp == 0      2.825e-02  2.601e-02   1.086   0.9999   
# 1 - IFSp == 0       4.432e-02  2.378e-02   1.863   0.9419   
# 43 - IFSp == 0     -2.060e-03  2.672e-02  -0.077   1.0000   
# OP4 - IFSp == 0    -2.490e-02  2.464e-02  -1.011   1.0000   
# PFop - IFSp == 0   -1.431e-02  2.529e-02  -0.566   1.0000   
# PF - IFSp == 0     -4.611e-03  2.732e-02  -0.169   1.0000   
# PSL - IFSp == 0    -3.756e-02  3.196e-02  -1.175   0.9998   
# STV - IFSp == 0    -1.415e-02  3.025e-02  -0.468   1.0000   
# A5 - IFSp == 0     -3.937e-02  2.601e-02  -1.514   0.9936   
# A4 - IFSp == 0     -3.730e-03  2.433e-02  -0.153   1.0000   
# 6r - 44 == 0       -3.327e-02  2.524e-02  -1.318   0.9989   
# 6v - 44 == 0        1.740e-03  2.151e-02   0.081   1.0000   
# 55b - 44 == 0       9.357e-03  2.320e-02   0.403   1.0000   
# 4 - 44 == 0        -5.053e-03  2.524e-02  -0.200   1.0000   
# 3a - 44 == 0       -1.780e-02  2.652e-02  -0.671   1.0000   
# 3b - 44 == 0        2.696e-02  2.381e-02   1.132   0.9999   
# 1 - 44 == 0         4.303e-02  2.136e-02   2.014   0.8878   
# 43 - 44 == 0       -3.355e-03  2.459e-02  -0.136   1.0000   
# OP4 - 44 == 0      -2.619e-02  2.231e-02  -1.174   0.9998   
# PFop - 44 == 0     -1.560e-02  2.303e-02  -0.678   1.0000   
# PF - 44 == 0       -5.906e-03  2.524e-02  -0.234   1.0000   
# PSL - 44 == 0      -3.885e-02  3.019e-02  -1.287   0.9992   
# STV - 44 == 0      -1.544e-02  2.838e-02  -0.544   1.0000   
# A5 - 44 == 0       -4.067e-02  2.381e-02  -1.708   0.9749   
# A4 - 44 == 0       -5.026e-03  2.197e-02  -0.229   1.0000   
# 6v - 6r == 0        3.501e-02  2.106e-02   1.662   0.9812   
# 55b - 6r == 0       4.263e-02  2.278e-02   1.871   0.9399   
# 4 - 6r == 0         2.822e-02  2.486e-02   1.135   0.9999   
# 3a - 6r == 0        1.547e-02  2.615e-02   0.591   1.0000   
# 3b - 6r == 0        6.023e-02  2.340e-02   2.574   0.5135   
# 1 - 6r == 0         7.630e-02  2.090e-02   3.650   0.0371 * 
#   43 - 6r == 0        2.992e-02  2.419e-02   1.237   0.9995   
# OP4 - 6r == 0       7.077e-03  2.187e-02   0.324   1.0000   
# PFop - 6r == 0      1.767e-02  2.260e-02   0.782   1.0000   
# PF - 6r == 0        2.736e-02  2.486e-02   1.101   0.9999   
# PSL - 6r == 0      -5.579e-03  2.987e-02  -0.187   1.0000   
# STV - 6r == 0       1.783e-02  2.804e-02   0.636   1.0000   
# A5 - 6r == 0       -7.396e-03  2.340e-02  -0.316   1.0000   
# A4 - 6r == 0        2.825e-02  2.153e-02   1.312   0.9989   
# 55b - 6v == 0       7.617e-03  1.857e-02   0.410   1.0000   
# 4 - 6v == 0        -6.792e-03  2.106e-02  -0.323   1.0000   
# 3a - 6v == 0       -1.954e-02  2.258e-02  -0.866   1.0000   
# 3b - 6v == 0        2.522e-02  1.932e-02   1.305   0.9990   
# 1 - 6v == 0         4.129e-02  1.621e-02   2.547   0.5334   
# 43 - 6v == 0       -5.095e-03  2.027e-02  -0.251   1.0000   
# OP4 - 6v == 0      -2.793e-02  1.744e-02  -1.602   0.9876   
# PFop - 6v == 0     -1.734e-02  1.835e-02  -0.945   1.0000   
# PF - 6v == 0       -7.646e-03  2.106e-02  -0.363   1.0000   
# PSL - 6v == 0      -4.059e-02  2.680e-02  -1.515   0.9935   
# STV - 6v == 0      -1.718e-02  2.474e-02  -0.694   1.0000   
# A5 - 6v == 0       -4.241e-02  1.932e-02  -2.195   0.7913   
# A4 - 6v == 0       -6.765e-03  1.700e-02  -0.398   1.0000   
# 4 - 55b == 0       -1.441e-02  2.278e-02  -0.633   1.0000   
# 3a - 55b == 0      -2.716e-02  2.419e-02  -1.123   0.9999   
# 3b - 55b == 0       1.760e-02  2.118e-02   0.831   1.0000   
# 1 - 55b == 0        3.367e-02  1.839e-02   1.831   0.9507   
# 43 - 55b == 0      -1.271e-02  2.206e-02  -0.576   1.0000   
# OP4 - 55b == 0     -3.555e-02  1.948e-02  -1.825   0.9526   
# PFop - 55b == 0    -2.496e-02  2.030e-02  -1.230   0.9996   
# PF - 55b == 0      -1.526e-02  2.278e-02  -0.670   1.0000   
# PSL - 55b == 0     -4.821e-02  2.817e-02  -1.711   0.9747   
# STV - 55b == 0     -2.480e-02  2.622e-02  -0.946   1.0000   
# A5 - 55b == 0      -5.002e-02  2.118e-02  -2.361   0.6755   
# A4 - 55b == 0      -1.438e-02  1.909e-02  -0.753   1.0000   
# 3a - 4 == 0        -1.275e-02  2.615e-02  -0.488   1.0000   
# 3b - 4 == 0         3.201e-02  2.340e-02   1.368   0.9982   
# 1 - 4 == 0          4.808e-02  2.090e-02   2.300   0.7217   
# 43 - 4 == 0         1.698e-03  2.419e-02   0.070   1.0000   
# OP4 - 4 == 0       -2.114e-02  2.187e-02  -0.967   1.0000   
# PFop - 4 == 0      -1.055e-02  2.260e-02  -0.467   1.0000   
# PF - 4 == 0        -8.534e-04  2.486e-02  -0.034   1.0000   
# PSL - 4 == 0       -3.380e-02  2.987e-02  -1.131   0.9999   
# STV - 4 == 0       -1.039e-02  2.804e-02  -0.370   1.0000   
# A5 - 4 == 0        -3.561e-02  2.340e-02  -1.522   0.9931   
# A4 - 4 == 0         2.708e-05  2.153e-02   0.001   1.0000   
# 3b - 3a == 0        4.476e-02  2.478e-02   1.807   0.9570   
# 1 - 3a == 0         6.083e-02  2.243e-02   2.712   0.4089   
# 43 - 3a == 0        1.445e-02  2.552e-02   0.566   1.0000   
# OP4 - 3a == 0      -8.390e-03  2.333e-02  -0.360   1.0000   
# PFop - 3a == 0      2.203e-03  2.402e-02   0.092   1.0000   
# PF - 3a == 0        1.190e-02  2.615e-02   0.455   1.0000   
# PSL - 3a == 0      -2.105e-02  3.096e-02  -0.680   1.0000   
# STV - 3a == 0       2.363e-03  2.920e-02   0.081   1.0000   
# A5 - 3a == 0       -2.286e-02  2.478e-02  -0.923   1.0000   
# A4 - 3a == 0        1.278e-02  2.301e-02   0.555   1.0000   
# 1 - 3b == 0         1.607e-02  1.915e-02   0.839   1.0000   
# 43 - 3b == 0       -3.031e-02  2.270e-02  -1.336   0.9987   
# OP4 - 3b == 0      -5.315e-02  2.020e-02  -2.631   0.4674   
# PFop - 3b == 0     -4.256e-02  2.099e-02  -2.027   0.8823   
# PF - 3b == 0       -3.287e-02  2.340e-02  -1.404   0.9974   
# PSL - 3b == 0      -6.581e-02  2.867e-02  -2.295   0.7232   
# STV - 3b == 0      -4.240e-02  2.676e-02  -1.584   0.9888   
# A5 - 3b == 0       -6.763e-02  2.185e-02  -3.095   0.1799   
# A4 - 3b == 0       -3.199e-02  1.983e-02  -1.613   0.9867   
# 43 - 1 == 0        -4.638e-02  2.011e-02  -2.306   0.7161   
# OP4 - 1 == 0       -6.922e-02  1.725e-02  -4.013   0.0100 * 
#   PFop - 1 == 0      -5.863e-02  1.817e-02  -3.227   0.1276   
# PF - 1 == 0        -4.893e-02  2.090e-02  -2.341   0.6927   
# PSL - 1 == 0       -8.188e-02  2.668e-02  -3.069   0.1908   
# STV - 1 == 0       -5.847e-02  2.461e-02  -2.376   0.6647   
# A5 - 1 == 0        -8.369e-02  1.915e-02  -4.370    <0.01 **
#   A4 - 1 == 0        -4.805e-02  1.681e-02  -2.859   0.3087   
# OP4 - 43 == 0      -2.284e-02  2.111e-02  -1.082   0.9999   
# PFop - 43 == 0     -1.225e-02  2.187e-02  -0.560   1.0000   
# PF - 43 == 0       -2.551e-03  2.419e-02  -0.105   1.0000   
# PSL - 43 == 0      -3.550e-02  2.932e-02  -1.210   0.9997   
# STV - 43 == 0      -1.209e-02  2.746e-02  -0.440   1.0000   
# A5 - 43 == 0       -3.731e-02  2.270e-02  -1.644   0.9835   
# A4 - 43 == 0       -1.671e-03  2.076e-02  -0.080   1.0000   
# PFop - OP4 == 0     1.059e-02  1.927e-02   0.550   1.0000   
# PF - OP4 == 0       2.029e-02  2.187e-02   0.928   1.0000   
# PSL - OP4 == 0     -1.266e-02  2.744e-02  -0.461   1.0000   
# STV - OP4 == 0      1.075e-02  2.543e-02   0.423   1.0000   
# A5 - OP4 == 0      -1.447e-02  2.020e-02  -0.716   1.0000   
# A4 - OP4 == 0       2.117e-02  1.800e-02   1.176   0.9998   
# PF - PFop == 0      9.695e-03  2.260e-02   0.429   1.0000   
# PSL - PFop == 0    -2.325e-02  2.803e-02  -0.830   1.0000   
# STV - PFop == 0     1.602e-04  2.607e-02   0.006   1.0000   
# A5 - PFop == 0     -2.507e-02  2.099e-02  -1.194   0.9997   
# A4 - PFop == 0      1.058e-02  1.888e-02   0.560   1.0000   
# PSL - PF == 0      -3.294e-02  2.987e-02  -1.103   0.9999   
# STV - PF == 0      -9.535e-03  2.804e-02  -0.340   1.0000   
# A5 - PF == 0       -3.476e-02  2.340e-02  -1.485   0.9948   
# A4 - PF == 0        8.805e-04  2.153e-02   0.041   1.0000   
# STV - PSL == 0      2.341e-02  3.257e-02   0.719   1.0000   
# A5 - PSL == 0      -1.816e-03  2.867e-02  -0.063   1.0000   
# A4 - PSL == 0       3.382e-02  2.717e-02   1.245   0.9995   
# A5 - STV == 0      -2.523e-02  2.676e-02  -0.943   1.0000   
# A4 - STV == 0       1.042e-02  2.514e-02   0.414   1.0000   
# A4 - A5 == 0        3.564e-02  1.983e-02   1.797   0.9586   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)


aper_ecog_roi_agg %>%
  group_by(region,region1) %>%
  summarise(mean_exponent_s = mean(exponent_s,na.rm=TRUE), se_exponent_s=sd(exponent_s,na.rm=TRUE)/sqrt(n())) %>%
  ggplot() + 
  aes(x=region,y=mean_exponent_s,color=region1,ymax=mean_exponent_s+se_exponent_s,ymin=mean_exponent_s-se_exponent_s)+
  geom_errorbar() +
  geom_point()+
  scale_color_brewer(type='qual',palette='Set2') +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  theme(legend.position='top',legend.title = element_blank(),
        axis.title.x=element_blank(),legend.text=element_text(size=6)) +
  labs(y="Exponent\n region effect")
ggsave('fig/B12_28c_aper_exponent_vs_anatomy.pdf',width=5,height=4)

amod <- aov(exponent_s ~ region, data=aper_ecog_roi_agg)
amod
summary(amod)
tuk <- glht(amod, linfct = mcp(region = "Tukey"))
summary(tuk) 

# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: aov(formula = exponent_s ~ region, data = aper_ecog_roi_agg)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
# 8C - p9-46v == 0    0.019887   0.099764   0.199   1.0000    
# 8Av - p9-46v == 0  -0.104658   0.102006  -1.026   1.0000    
# IFSp - p9-46v == 0 -0.076266   0.112772  -0.676   1.0000    
# 44 - p9-46v == 0   -0.199575   0.106171  -1.880   0.9373    
# 6r - p9-46v == 0   -0.183548   0.104963  -1.749   0.9687    
# 6v - p9-46v == 0    0.066172   0.093357   0.709   1.0000    
# 55b - p9-46v == 0   0.076087   0.098547   0.772   1.0000    
# 4 - p9-46v == 0    -0.031352   0.104963  -0.299   1.0000    
# 3a - p9-46v == 0   -0.182024   0.109049  -1.669   0.9803    
# 3b - p9-46v == 0    0.053971   0.100448   0.537   1.0000    
# 1 - p9-46v == 0     0.015318   0.092891   0.165   1.0000    
# 43 - p9-46v == 0   -0.115018   0.102897  -1.118   0.9999    
# OP4 - p9-46v == 0  -0.234242   0.095783  -2.446   0.6122    
# PFop - p9-46v == 0 -0.154755   0.098004  -1.579   0.9895    
# PF - p9-46v == 0   -0.153111   0.104963  -1.459   0.9959    
# PSL - p9-46v == 0  -0.279195   0.121014  -2.307   0.7158    
# STV - p9-46v == 0  -0.166561   0.115078  -1.447   0.9962    
# A5 - p9-46v == 0   -0.162355   0.100448  -1.616   0.9863    
# A4 - p9-46v == 0   -0.176873   0.094747  -1.867   0.9404    
# 8Av - 8C == 0      -0.124545   0.080526  -1.547   0.9918    
# IFSp - 8C == 0     -0.096153   0.093790  -1.025   1.0000    
# 44 - 8C == 0       -0.219462   0.085741  -2.560   0.5241    
# 6r - 8C == 0       -0.203435   0.084241  -2.415   0.6351    
# 6v - 8C == 0        0.046285   0.069243   0.668   1.0000    
# 55b - 8C == 0       0.056200   0.076097   0.739   1.0000    
# 4 - 8C == 0        -0.051239   0.084241  -0.608   1.0000    
# 3a - 8C == 0       -0.201911   0.089280  -2.262   0.7473    
# 3b - 8C == 0        0.034083   0.078544   0.434   1.0000    
# 1 - 8C == 0        -0.004569   0.068614  -0.067   1.0000    
# 43 - 8C == 0       -0.134905   0.081651  -1.652   0.9825    
# OP4 - 8C == 0      -0.254129   0.072482  -3.506   0.0567 .  
# PFop - 8C == 0     -0.174642   0.075392  -2.316   0.7083    
# PF - 8C == 0       -0.172998   0.084241  -2.054   0.8697    
# PSL - 8C == 0      -0.299082   0.103555  -2.888   0.2887    
# STV - 8C == 0      -0.186448   0.096551  -1.931   0.9214    
# A5 - 8C == 0       -0.182243   0.078544  -2.320   0.7065    
# A4 - 8C == 0       -0.196760   0.071107  -2.767   0.3660    
# IFSp - 8Av == 0     0.028392   0.096172   0.295   1.0000    
# 44 - 8Av == 0      -0.094917   0.088340  -1.074   0.9999    
# 6r - 8Av == 0      -0.078890   0.086884  -0.908   1.0000    
# 6v - 8Av == 0       0.170830   0.072437   2.358   0.6783    
# 55b - 8Av == 0      0.180745   0.079013   2.288   0.7287    
# 4 - 8Av == 0        0.073306   0.086884   0.844   1.0000    
# 3a - 8Av == 0      -0.077366   0.091778  -0.843   1.0000    
# 3b - 8Av == 0       0.158629   0.081372   1.949   0.9144    
# 1 - 8Av == 0        0.119976   0.071835   1.670   0.9806    
# 43 - 8Av == 0      -0.010360   0.084376  -0.123   1.0000    
# OP4 - 8Av == 0     -0.129584   0.075539  -1.715   0.9740    
# PFop - 8Av == 0    -0.050097   0.078335  -0.640   1.0000    
# PF - 8Av == 0      -0.048453   0.086884  -0.558   1.0000    
# PSL - 8Av == 0     -0.174537   0.105716  -1.651   0.9826    
# STV - 8Av == 0     -0.061903   0.098866  -0.626   1.0000    
# A5 - 8Av == 0      -0.057697   0.081372  -0.709   1.0000    
# A4 - 8Av == 0      -0.072215   0.074220  -0.973   1.0000    
# 44 - IFSp == 0     -0.123309   0.100579  -1.226   0.9996    
# 6r - IFSp == 0     -0.107282   0.099303  -1.080   0.9999    
# 6v - IFSp == 0      0.142438   0.086944   1.638   0.9841    
# 55b - IFSp == 0     0.152353   0.092495   1.647   0.9832    
# 4 - IFSp == 0       0.044914   0.099303   0.452   1.0000    
# 3a - IFSp == 0     -0.105758   0.103612  -1.021   1.0000    
# 3b - IFSp == 0      0.130236   0.094519   1.378   0.9980    
# 1 - IFSp == 0       0.091584   0.086444   1.059   1.0000    
# 43 - IFSp == 0     -0.038752   0.097116  -0.399   1.0000    
# OP4 - IFSp == 0    -0.157976   0.089545  -1.764   0.9656    
# PFop - IFSp == 0   -0.078490   0.091917  -0.854   1.0000    
# PF - IFSp == 0     -0.076845   0.099303  -0.774   1.0000    
# PSL - IFSp == 0    -0.202930   0.116139  -1.747   0.9688    
# STV - IFSp == 0    -0.090295   0.109940  -0.821   1.0000    
# A5 - IFSp == 0     -0.086090   0.094519  -0.911   1.0000    
# A4 - IFSp == 0     -0.100607   0.088436  -1.138   0.9999    
# 6r - 44 == 0        0.016027   0.091739   0.175   1.0000    
# 6v - 44 == 0        0.265748   0.078193   3.399   0.0784 .  
# 55b - 44 == 0       0.275662   0.084322   3.269   0.1145    
# 4 - 44 == 0         0.168223   0.091739   1.834   0.9500    
# 3a - 44 == 0        0.017551   0.096386   0.182   1.0000    
# 3b - 44 == 0        0.253546   0.086537   2.930   0.2647    
# 1 - 44 == 0         0.214893   0.077636   2.768   0.3682    
# 43 - 44 == 0        0.084557   0.089367   0.946   1.0000    
# OP4 - 44 == 0      -0.034667   0.081075  -0.428   1.0000    
# PFop - 44 == 0      0.044820   0.083687   0.536   1.0000    
# PF - 44 == 0        0.046464   0.091739   0.506   1.0000    
# PSL - 44 == 0      -0.079620   0.109741  -0.726   1.0000    
# STV - 44 == 0       0.033015   0.103158   0.320   1.0000    
# A5 - 44 == 0        0.037220   0.086537   0.430   1.0000    
# A4 - 44 == 0        0.022702   0.079848   0.284   1.0000    
# 6v - 6r == 0        0.249721   0.076545   3.262   0.1143    
# 55b - 6r == 0       0.259635   0.082796   3.136   0.1619    
# 4 - 6r == 0         0.152196   0.090338   1.685   0.9785    
# 3a - 6r == 0        0.001524   0.095054   0.016   1.0000    
# 3b - 6r == 0        0.237519   0.085050   2.793   0.3512    
# 1 - 6r == 0         0.198866   0.075976   2.617   0.4769    
# 43 - 6r == 0        0.068530   0.087928   0.779   1.0000    
# OP4 - 6r == 0      -0.050694   0.079487  -0.638   1.0000    
# PFop - 6r == 0      0.028793   0.082149   0.350   1.0000    
# PF - 6r == 0        0.030437   0.090338   0.337   1.0000    
# PSL - 6r == 0      -0.095647   0.108573  -0.881   1.0000    
# STV - 6r == 0       0.016988   0.101915   0.167   1.0000    
# A5 - 6r == 0        0.021193   0.085050   0.249   1.0000    
# A4 - 6r == 0        0.006675   0.078235   0.085   1.0000    
# 55b - 6v == 0       0.009915   0.067479   0.147   1.0000    
# 4 - 6v == 0        -0.097525   0.076545  -1.274   0.9993    
# 3a - 6v == 0       -0.248196   0.082058  -3.025   0.2133    
# 3b - 6v == 0       -0.012202   0.070226  -0.174   1.0000    
# 1 - 6v == 0        -0.050854   0.058912  -0.863   1.0000    
# 43 - 6v == 0       -0.181190   0.073686  -2.459   0.6017    
# OP4 - 6v == 0      -0.300415   0.063375  -4.740    <0.01 ***
#   PFop - 6v == 0     -0.220928   0.066683  -3.313   0.1007    
# PF - 6v == 0       -0.219283   0.076545  -2.865   0.3055    
# PSL - 6v == 0      -0.345368   0.097397  -3.546   0.0508 .  
# STV - 6v == 0      -0.232733   0.089916  -2.588   0.4997    
# A5 - 6v == 0       -0.228528   0.070226  -3.254   0.1190    
# A4 - 6v == 0       -0.243045   0.061797  -3.933   0.0131 *  
#   4 - 55b == 0       -0.107439   0.082796  -1.298   0.9991    
# 3a - 55b == 0      -0.258111   0.087918  -2.936   0.2602    
# 3b - 55b == 0      -0.022116   0.076992  -0.287   1.0000    
# 1 - 55b == 0       -0.060769   0.066833  -0.909   1.0000    
# 43 - 55b == 0      -0.191105   0.080160  -2.384   0.6591    
# OP4 - 55b == 0     -0.310329   0.070798  -4.383    <0.01 ** 
#   PFop - 55b == 0    -0.230842   0.073775  -3.129   0.1656    
# PF - 55b == 0      -0.229198   0.082796  -2.768   0.3670    
# PSL - 55b == 0     -0.355282   0.102383  -3.470   0.0650 .  
# STV - 55b == 0     -0.242648   0.095294  -2.546   0.5345    
# A5 - 55b == 0      -0.238442   0.076992  -3.097   0.1798    
# A4 - 55b == 0      -0.252960   0.069390  -3.645   0.0373 *  
#   3a - 4 == 0        -0.150672   0.095054  -1.585   0.9889    
# 3b - 4 == 0         0.085323   0.085050   1.003   1.0000    
# 1 - 4 == 0          0.046670   0.075976   0.614   1.0000    
# 43 - 4 == 0        -0.083666   0.087928  -0.952   1.0000    
# OP4 - 4 == 0       -0.202890   0.079487  -2.553   0.5279    
# PFop - 4 == 0      -0.123403   0.082149  -1.502   0.9941    
# PF - 4 == 0        -0.121759   0.090338  -1.348   0.9985    
# PSL - 4 == 0       -0.247843   0.108573  -2.283   0.7317    
# STV - 4 == 0       -0.135208   0.101915  -1.327   0.9988    
# A5 - 4 == 0        -0.131003   0.085050  -1.540   0.9920    
# A4 - 4 == 0        -0.145521   0.078235  -1.860   0.9433    
# 3b - 3a == 0        0.235995   0.090044   2.621   0.4774    
# 1 - 3a == 0         0.197342   0.081528   2.421   0.6321    
# 43 - 3a == 0        0.067006   0.092767   0.722   1.0000    
# OP4 - 3a == 0      -0.052218   0.084809  -0.616   1.0000    
# PFop - 3a == 0      0.027269   0.087309   0.312   1.0000    
# PF - 3a == 0        0.028913   0.095054   0.304   1.0000    
# PSL - 3a == 0      -0.097171   0.112527  -0.864   1.0000    
# STV - 3a == 0       0.015463   0.106118   0.146   1.0000    
# A5 - 3a == 0        0.019669   0.090044   0.218   1.0000    
# A4 - 3a == 0        0.005151   0.083637   0.062   1.0000    
# 1 - 3b == 0        -0.038652   0.069606  -0.555   1.0000    
# 43 - 3b == 0       -0.168989   0.082487  -2.049   0.8724    
# OP4 - 3b == 0      -0.288213   0.073422  -3.925   0.0130 *  
#   PFop - 3b == 0     -0.208726   0.076296  -2.736   0.3919    
# PF - 3b == 0       -0.207081   0.085050  -2.435   0.6191    
# PSL - 3b == 0      -0.333166   0.104214  -3.197   0.1379    
# STV - 3b == 0      -0.220531   0.097259  -2.267   0.7428    
# A5 - 3b == 0       -0.216326   0.079411  -2.724   0.3974    
# A4 - 3b == 0       -0.230844   0.072065  -3.203   0.1367    
# 43 - 1 == 0        -0.130336   0.073095  -1.783   0.9617    
# OP4 - 1 == 0       -0.249560   0.062686  -3.981   0.0117 *  
#   PFop - 1 == 0      -0.170073   0.066030  -2.576   0.5092    
# PF - 1 == 0        -0.168429   0.075976  -2.217   0.7768    
# PSL - 1 == 0       -0.294513   0.096951  -3.038   0.2058    
# STV - 1 == 0       -0.181879   0.089432  -2.034   0.8795    
# A5 - 1 == 0        -0.177674   0.069606  -2.553   0.5281    
# A4 - 1 == 0        -0.192191   0.061091  -3.146   0.1564    
# OP4 - 43 == 0      -0.119224   0.076737  -1.554   0.9912    
# PFop - 43 == 0     -0.039737   0.079492  -0.500   1.0000    
# PF - 43 == 0       -0.038093   0.087928  -0.433   1.0000    
# PSL - 43 == 0      -0.164177   0.106576  -1.540   0.9919    
# STV - 43 == 0      -0.051542   0.099785  -0.517   1.0000    
# A5 - 43 == 0       -0.047337   0.082487  -0.574   1.0000    
# A4 - 43 == 0       -0.061855   0.075440  -0.820   1.0000    
# PFop - OP4 == 0     0.079487   0.070041   1.135   0.9999    
# PF - OP4 == 0       0.081131   0.079487   1.021   1.0000    
# PSL - OP4 == 0     -0.044953   0.099726  -0.451   1.0000    
# STV - OP4 == 0      0.067682   0.092433   0.732   1.0000    
# A5 - OP4 == 0       0.071887   0.073422   0.979   1.0000    
# A4 - OP4 == 0       0.057369   0.065406   0.877   1.0000    
# PF - PFop == 0      0.001645   0.082149   0.020   1.0000    
# PSL - PFop == 0    -0.124440   0.101860  -1.222   0.9996    
# STV - PFop == 0    -0.011805   0.094732  -0.125   1.0000    
# A5 - PFop == 0     -0.007600   0.076296  -0.100   1.0000    
# A4 - PFop == 0     -0.022118   0.068617  -0.322   1.0000    
# PSL - PF == 0      -0.126085   0.108573  -1.161   0.9998    
# STV - PF == 0      -0.013450   0.101915  -0.132   1.0000    
# A5 - PF == 0       -0.009245   0.085050  -0.109   1.0000    
# A4 - PF == 0       -0.023762   0.078235  -0.304   1.0000    
# STV - PSL == 0      0.112635   0.118380   0.951   1.0000    
# A5 - PSL == 0       0.116840   0.104214   1.121   0.9999    
# A4 - PSL == 0       0.102322   0.098731   1.036   1.0000    
# A5 - STV == 0       0.004205   0.097259   0.043   1.0000    
# A4 - STV == 0      -0.010313   0.091359  -0.113   1.0000    
# A4 - A5 == 0       -0.014518   0.072065  -0.201   1.0000    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)




aper_ecog_roi_agg %>%
  group_by(region,region1) %>%
  summarise(mean_offset_s = mean(offset_s,na.rm=TRUE), se_offset_s=sd(offset_s,na.rm=TRUE)/sqrt(n())) %>%
  ggplot() + 
  aes(x=region,y=mean_offset_s,color=region1,ymax=mean_offset_s+se_offset_s,ymin=mean_offset_s-se_offset_s)+
  geom_errorbar() +
  geom_point()+
  scale_color_brewer(type='qual',palette='Set2') +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  theme(legend.position='top',legend.title = element_blank(),
        axis.title.x=element_blank(),legend.text=element_text(size=6)) +
  labs(y="Offset\n region effect")
ggsave('fig/B12_28c_aper_offset_vs_anatomy.pdf',width=5,height=4)

amod <- aov(offset_s ~ region, data=aper_ecog_roi_agg)
amod
summary(amod)
tuk <- glht(amod, linfct = mcp(region = "Tukey"))
summary(tuk) 

# 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: aov(formula = offset_s ~ region, data = aper_ecog_roi_agg)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
# 8C - p9-46v == 0   -0.090979   0.106538  -0.854   1.0000    
# 8Av - p9-46v == 0  -0.053797   0.108933  -0.494   1.0000    
# IFSp - p9-46v == 0 -0.165943   0.120430  -1.378   0.9980    
# 44 - p9-46v == 0   -0.413911   0.113381  -3.651   0.0354 *  
#   6r - p9-46v == 0   -0.393780   0.112091  -3.513   0.0554 .  
# 6v - p9-46v == 0   -0.074302   0.099696  -0.745   1.0000    
# 55b - p9-46v == 0  -0.038524   0.105239  -0.366   1.0000    
# 4 - p9-46v == 0    -0.152022   0.112091  -1.356   0.9984    
# 3a - p9-46v == 0   -0.471299   0.116454  -4.047    <0.01 ** 
#   3b - p9-46v == 0   -0.149297   0.107270  -1.392   0.9977    
# 1 - p9-46v == 0    -0.251768   0.099199  -2.538   0.5392    
# 43 - p9-46v == 0   -0.235702   0.109884  -2.145   0.8224    
# OP4 - p9-46v == 0  -0.421441   0.102288  -4.120    <0.01 ** 
#   PFop - p9-46v == 0 -0.259350   0.104659  -2.478   0.5866    
# PF - p9-46v == 0   -0.250946   0.112091  -2.239   0.7628    
# PSL - p9-46v == 0  -0.207374   0.129231  -1.605   0.9872    
# STV - p9-46v == 0  -0.212652   0.122893  -1.730   0.9716    
# A5 - p9-46v == 0   -0.177860   0.107270  -1.658   0.9819    
# A4 - p9-46v == 0   -0.282794   0.101181  -2.795   0.3503    
# 8Av - 8C == 0       0.037183   0.085994   0.432   1.0000    
# IFSp - 8C == 0     -0.074964   0.100159  -0.748   1.0000    
# 44 - 8C == 0       -0.322932   0.091563  -3.527   0.0543 .  
# 6r - 8C == 0       -0.302801   0.089961  -3.366   0.0888 .  
# 6v - 8C == 0        0.016678   0.073946   0.226   1.0000    
# 55b - 8C == 0       0.052456   0.081264   0.645   1.0000    
# 4 - 8C == 0        -0.061043   0.089961  -0.679   1.0000    
# 3a - 8C == 0       -0.380319   0.095342  -3.989   0.0103 *  
#   3b - 8C == 0       -0.058318   0.083877  -0.695   1.0000    
# 1 - 8C == 0        -0.160789   0.073274  -2.194   0.7910    
# 43 - 8C == 0       -0.144722   0.087196  -1.660   0.9816    
# OP4 - 8C == 0      -0.330461   0.077404  -4.269    <0.01 ** 
#   PFop - 8C == 0     -0.168371   0.080512  -2.091   0.8509    
# PF - 8C == 0       -0.159967   0.089961  -1.778   0.9628    
# PSL - 8C == 0      -0.116395   0.110587  -1.053   1.0000    
# STV - 8C == 0      -0.121673   0.103108  -1.180   0.9998    
# A5 - 8C == 0       -0.086881   0.083877  -1.036   1.0000    
# A4 - 8C == 0       -0.191814   0.075936  -2.526   0.5503    
# IFSp - 8Av == 0    -0.112146   0.102703  -1.092   0.9999    
# 44 - 8Av == 0      -0.360114   0.094339  -3.817   0.0193 *  
#   6r - 8Av == 0      -0.339983   0.092784  -3.664   0.0347 *  
#   6v - 8Av == 0      -0.020505   0.077355  -0.265   1.0000    
# 55b - 8Av == 0      0.015273   0.084379   0.181   1.0000    
# 4 - 8Av == 0       -0.098225   0.092784  -1.059   1.0000    
# 3a - 8Av == 0      -0.417502   0.098011  -4.260    <0.01 ** 
#   3b - 8Av == 0      -0.095500   0.086898  -1.099   0.9999    
# 1 - 8Av == 0       -0.197971   0.076713  -2.581   0.5070    
# 43 - 8Av == 0      -0.181905   0.090106  -2.019   0.8868    
# OP4 - 8Av == 0     -0.367644   0.080668  -4.557    <0.01 ***
#   PFop - 8Av == 0    -0.205553   0.083655  -2.457   0.6033    
# PF - 8Av == 0      -0.197150   0.092784  -2.125   0.8320    
# PSL - 8Av == 0     -0.153577   0.112895  -1.360   0.9983    
# STV - 8Av == 0     -0.158855   0.105580  -1.505   0.9939    
# A5 - 8Av == 0      -0.124063   0.086898  -1.428   0.9968    
# A4 - 8Av == 0      -0.228997   0.079260  -2.889   0.2886    
# 44 - IFSp == 0     -0.247968   0.107409  -2.309   0.7145    
# 6r - IFSp == 0     -0.227837   0.106047  -2.148   0.8190    
# 6v - IFSp == 0      0.091642   0.092849   0.987   1.0000    
# 55b - IFSp == 0     0.127419   0.098776   1.290   0.9992    
# 4 - IFSp == 0       0.013921   0.106047   0.131   1.0000    
# 3a - IFSp == 0     -0.305355   0.110648  -2.760   0.3734    
# 3b - IFSp == 0      0.016646   0.100937   0.165   1.0000    
# 1 - IFSp == 0      -0.085825   0.092314  -0.930   1.0000    
# 43 - IFSp == 0     -0.069758   0.103711  -0.673   1.0000    
# OP4 - IFSp == 0    -0.255497   0.095626  -2.672   0.4362    
# PFop - IFSp == 0   -0.093407   0.098158  -0.952   1.0000    
# PF - IFSp == 0     -0.085003   0.106047  -0.802   1.0000    
# PSL - IFSp == 0    -0.041431   0.124025  -0.334   1.0000    
# STV - IFSp == 0    -0.046709   0.117406  -0.398   1.0000    
# A5 - IFSp == 0     -0.011917   0.100937  -0.118   1.0000    
# A4 - IFSp == 0     -0.116850   0.094441  -1.237   0.9995    
# 6r - 44 == 0        0.020131   0.097968   0.205   1.0000    
# 6v - 44 == 0        0.339609   0.083503   4.067    <0.01 ** 
#   55b - 44 == 0       0.375387   0.090048   4.169    <0.01 ** 
#   4 - 44 == 0         0.261889   0.097968   2.673   0.4369    
# 3a - 44 == 0       -0.057388   0.102932  -0.558   1.0000    
# 3b - 44 == 0        0.264614   0.092413   2.863   0.3043    
# 1 - 44 == 0         0.162143   0.082908   1.956   0.9125    
# 43 - 44 == 0        0.178210   0.095436   1.867   0.9412    
# OP4 - 44 == 0      -0.007530   0.086581  -0.087   1.0000    
# PFop - 44 == 0      0.154561   0.089370   1.729   0.9721    
# PF - 44 == 0        0.162965   0.097968   1.663   0.9812    
# PSL - 44 == 0       0.206537   0.117193   1.762   0.9660    
# STV - 44 == 0       0.201259   0.110164   1.827   0.9516    
# A5 - 44 == 0        0.236051   0.092413   2.554   0.5262    
# A4 - 44 == 0        0.131117   0.085270   1.538   0.9921    
# 6v - 6r == 0        0.319478   0.081743   3.908   0.0146 *  
#   55b - 6r == 0       0.355256   0.088418   4.018    <0.01 ** 
#   4 - 6r == 0         0.241758   0.096472   2.506   0.5637    
# 3a - 6r == 0       -0.077519   0.101509  -0.764   1.0000    
# 3b - 6r == 0        0.244483   0.090826   2.692   0.4239    
# 1 - 6r == 0         0.142012   0.081135   1.750   0.9685    
# 43 - 6r == 0        0.158079   0.093899   1.683   0.9786    
# OP4 - 6r == 0      -0.027661   0.084884  -0.326   1.0000    
# PFop - 6r == 0      0.134430   0.087727   1.532   0.9924    
# PF - 6r == 0        0.142834   0.096472   1.481   0.9950    
# PSL - 6r == 0       0.186406   0.115945   1.608   0.9869    
# STV - 6r == 0       0.181128   0.108836   1.664   0.9813    
# A5 - 6r == 0        0.215920   0.090826   2.377   0.6636    
# A4 - 6r == 0        0.110986   0.083548   1.328   0.9987    
# 55b - 6v == 0       0.035778   0.072061   0.496   1.0000    
# 4 - 6v == 0        -0.077720   0.081743  -0.951   1.0000    
# 3a - 6v == 0       -0.396997   0.087630  -4.530    <0.01 ** 
#   3b - 6v == 0       -0.074995   0.074995  -1.000   1.0000    
# 1 - 6v == 0        -0.177467   0.062912  -2.821   0.3310    
# 43 - 6v == 0       -0.161400   0.078690  -2.051   0.8713    
# OP4 - 6v == 0      -0.347139   0.067678  -5.129    <0.01 ***
#   PFop - 6v == 0     -0.185048   0.071211  -2.599   0.4927    
# PF - 6v == 0       -0.176645   0.081743  -2.161   0.8119    
# PSL - 6v == 0      -0.133072   0.104011  -1.279   0.9992    
# STV - 6v == 0      -0.138350   0.096022  -1.441   0.9964    
# A5 - 6v == 0       -0.103558   0.074995  -1.381   0.9979    
# A4 - 6v == 0       -0.208492   0.065994  -3.159   0.1532    
# 4 - 55b == 0       -0.113498   0.088418  -1.284   0.9992    
# 3a - 55b == 0      -0.432775   0.093888  -4.609    <0.01 ***
#   3b - 55b == 0      -0.110773   0.082221  -1.347   0.9985    
# 1 - 55b == 0       -0.213244   0.071371  -2.988   0.2321    
# 43 - 55b == 0      -0.197178   0.085604  -2.303   0.7180    
# OP4 - 55b == 0     -0.382917   0.075606  -5.065    <0.01 ***
#   PFop - 55b == 0    -0.220826   0.078785  -2.803   0.3446    
# PF - 55b == 0      -0.212422   0.088418  -2.402   0.6454    
# PSL - 55b == 0     -0.168850   0.109335  -1.544   0.9917    
# STV - 55b == 0     -0.174128   0.101765  -1.711   0.9747    
# A5 - 55b == 0      -0.139336   0.082221  -1.695   0.9768    
# A4 - 55b == 0      -0.244270   0.074102  -3.296   0.1086    
# 3a - 4 == 0        -0.319277   0.101509  -3.145   0.1581    
# 3b - 4 == 0         0.002725   0.090826   0.030   1.0000    
# 1 - 4 == 0         -0.099746   0.081135  -1.229   0.9996    
# 43 - 4 == 0        -0.083680   0.093899  -0.891   1.0000    
# OP4 - 4 == 0       -0.269419   0.084884  -3.174   0.1481    
# PFop - 4 == 0      -0.107328   0.087727  -1.223   0.9996    
# PF - 4 == 0        -0.098924   0.096472  -1.025   1.0000    
# PSL - 4 == 0       -0.055352   0.115945  -0.477   1.0000    
# STV - 4 == 0       -0.060630   0.108836  -0.557   1.0000    
# A5 - 4 == 0        -0.025838   0.090826  -0.284   1.0000    
# A4 - 4 == 0        -0.130772   0.083548  -1.565   0.9904    
# 3b - 3a == 0        0.322002   0.096159   3.349   0.0927 .  
# 1 - 3a == 0         0.219530   0.087064   2.521   0.5531    
# 43 - 3a == 0        0.235597   0.099067   2.378   0.6624    
# OP4 - 3a == 0       0.049858   0.090568   0.551   1.0000    
# PFop - 3a == 0      0.211949   0.093238   2.273   0.7399    
# PF - 3a == 0        0.220352   0.101509   2.171   0.8066    
# PSL - 3a == 0       0.263925   0.120169   2.196   0.7903    
# STV - 3a == 0       0.258646   0.113324   2.282   0.7344    
# A5 - 3a == 0        0.293438   0.096159   3.052   0.2002    
# A4 - 3a == 0        0.188505   0.089316   2.111   0.8405    
# 1 - 3b == 0        -0.102471   0.074333  -1.379   0.9980    
# 43 - 3b == 0       -0.086404   0.088088  -0.981   1.0000    
# OP4 - 3b == 0      -0.272144   0.078408  -3.471   0.0641 .  
# PFop - 3b == 0     -0.110053   0.081477  -1.351   0.9984    
# PF - 3b == 0       -0.101649   0.090826  -1.119   0.9999    
# PSL - 3b == 0      -0.058077   0.111291  -0.522   1.0000    
# STV - 3b == 0      -0.063355   0.103863  -0.610   1.0000    
# A5 - 3b == 0       -0.028563   0.084804  -0.337   1.0000    
# A4 - 3b == 0       -0.133497   0.076958  -1.735   0.9711    
# 43 - 1 == 0         0.016067   0.078058   0.206   1.0000    
# OP4 - 1 == 0       -0.169672   0.066943  -2.535   0.5425    
# PFop - 1 == 0      -0.007582   0.070513  -0.108   1.0000    
# PF - 1 == 0         0.000822   0.081135   0.010   1.0000    
# PSL - 1 == 0        0.044394   0.103534   0.429   1.0000    
# STV - 1 == 0        0.039116   0.095505   0.410   1.0000    
# A5 - 1 == 0         0.073908   0.074333   0.994   1.0000    
# A4 - 1 == 0        -0.031026   0.065240  -0.476   1.0000    
# OP4 - 43 == 0      -0.185739   0.081948  -2.267   0.7441    
# PFop - 43 == 0     -0.023648   0.084890  -0.279   1.0000    
# PF - 43 == 0       -0.015245   0.093899  -0.162   1.0000    
# PSL - 43 == 0       0.028327   0.113813   0.249   1.0000    
# STV - 43 == 0       0.023049   0.106562   0.216   1.0000    
# A5 - 43 == 0        0.057841   0.088088   0.657   1.0000    
# A4 - 43 == 0       -0.047092   0.080563  -0.585   1.0000    
# PFop - OP4 == 0     0.162091   0.074797   2.167   0.8082    
# PF - OP4 == 0       0.170494   0.084884   2.009   0.8896    
# PSL - OP4 == 0      0.214067   0.106498   2.010   0.8893    
# STV - OP4 == 0      0.208788   0.098710   2.115   0.8379    
# A5 - OP4 == 0       0.243581   0.078408   3.107   0.1743    
# A4 - OP4 == 0       0.138647   0.069847   1.985   0.9000    
# PF - PFop == 0      0.008404   0.087727   0.096   1.0000    
# PSL - PFop == 0     0.051976   0.108777   0.478   1.0000    
# STV - PFop == 0     0.046698   0.101165   0.462   1.0000    
# A5 - PFop == 0      0.081490   0.081477   1.000   1.0000    
# A4 - PFop == 0     -0.023444   0.073276  -0.320   1.0000    
# PSL - PF == 0       0.043572   0.115945   0.376   1.0000    
# STV - PF == 0       0.038294   0.108836   0.352   1.0000    
# A5 - PF == 0        0.073086   0.090826   0.805   1.0000    
# A4 - PF == 0       -0.031848   0.083548  -0.381   1.0000    
# STV - PSL == 0     -0.005278   0.126418  -0.042   1.0000    
# A5 - PSL == 0       0.029514   0.111291   0.265   1.0000    
# A4 - PSL == 0      -0.075420   0.105435  -0.715   1.0000    
# A5 - STV == 0       0.034792   0.103863   0.335   1.0000    
# A4 - STV == 0      -0.070142   0.097563  -0.719   1.0000    
# A4 - A5 == 0       -0.104934   0.076958  -1.364   0.9982    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)


