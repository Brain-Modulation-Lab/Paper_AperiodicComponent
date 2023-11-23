library(coin)
library(multcomp)
library(lmerTest)
library(ggExtra)
library(ggradar)
library(tidyverse)


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

# DBS_aper <- read_csv('data/spectral_sensitivity_ap-params.csv') %>%
#   filter(freq_range=='1-200hz') %>%
#   select(-'...1') %>%
#   mutate(electrode=tolower(electrode)) %>%
#   left_join(subjects %>% dplyr::select(subject,dx,dbs_target)) %>%
#   mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
#   filter(!subsesel %in% artifacts$subsesel)
# 

elec_session <- read_tsv('data/electrode_session.txt') %>%
  dplyr::select(subject,session_id,electrode,type,side,
                mni_nonlinear_x,mni_nonlinear_y,mni_nonlinear_z,target,
                fs_anatomy,DISTAL_label_1,DISTAL_label_2,DISTAL_label_3,MOREL_label_1,MOREL_label_2,MOREL_label_3,HCPMMP1_label_1) %>%
  mutate(electrode=tolower(electrode)) 

DBS_aper <- read_csv('data/spectral_sensitivity_ap-params.csv') %>%
  filter(freq_range=='1-200hz') %>%
  select(-'...1') %>%
  select(-type) %>%
  left_join(elec_session) %>%
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

DBS_aper_coord_dx <- DBS_aper_filt %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel) %>%
  mutate(electrode_type = type, rsquared = R2)

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
  filter(rsquared > 0.95) %>%
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
  filter(n>=10)

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
# region1                        n p_off_total_fdr p_on_total_fdr p_off_total p_on_total r_off_total r_on_total
# <chr>                      <int>           <dbl>          <dbl>       <dbl>      <dbl>       <dbl>      <dbl>
#   1 auditory-associative        20           0.802          0.764       0.652     0.392      0.108        0.203
# 2 dorsolateral-prefrontal       14           0.854          0.242       0.748     0.0282     0.0946       0.584
# 3 inferior-frontal              10           0.764          0.242       0.498     0.0302     0.243        0.681
# 4 inferior-parietal             16           0.866          0.764       0.812     0.541     -0.0648       0.165
# 5 posterior-opercular           17           0.970          0.764       0.970     0.573      0.00981      0.147
# 6 premotor                      20           0.764          0.539       0.356     0.135      0.218        0.346
# 7 sensorimotor                  23           0.764          0.348       0.328     0.0652     0.214        0.391
# 8 temporo-parietal-occipital    12           0.764          0.751       0.514     0.235     -0.210       -0.371

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
  filter(n>=10)

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

# region1                     n p_on_axial_fdr p_on_bradykinesia_fdr p_on_rigidity_fdr p_on_tremor_fdr p_off_axial_fdr p_off_bradykinesia_fdr
# <chr>                   <int>          <dbl>                 <dbl>             <dbl>           <dbl>           <dbl>                  <dbl>
# 1 dorsolateral-prefrontal    14          0.714                0.0515             0.638           0.131           0.638                  0.182
# 2 inferior-frontal           10          0.728                0.0419             0.483           0.728           0.728                  0.154

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
ggsave('fig/B22_11a_exponent-vs-UPDRS-on-per-region.pdf',width=7,height=4.5)

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
ggsave('fig/B22_11b_exponent-vs-UPDRS-off-per-region.pdf',width=7,height=4.5)


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
ggsave('fig/B22_12a_exponent-vs-UPDRS-on-axial-per-region.pdf',width=4,height=3)

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
ggsave('fig/B22_12a_exponent-vs-UPDRS-on-bradykinesia-per-region.pdf',width=4,height=3)

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
ggsave('fig/B22_12a_exponent-vs-UPDRS-on-rigidity-per-region.pdf',width=4,height=3)

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
ggsave('fig/B22_12a_exponent-vs-UPDRS-on-tremor-per-region.pdf',width=4,height=3)


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
ggsave('fig/B22_12a_exponent-vs-UPDRS-off-axial-per-region.pdf',width=4,height=3)

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
ggsave('fig/B22_12a_exponent-vs-UPDRS-off-bradykinesia-per-region.pdf',width=4,height=3)

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
ggsave('fig/B22_12a_exponent-vs-UPDRS-off-rigidity-per-region.pdf',width=4,height=3)

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
ggsave('fig/B22_12a_exponent-vs-UPDRS-off-tremor-per-region.pdf',width=4,height=3)



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
  
  ggsave(paste0('fig/B22_13a_radar-plot_',ur,'.pdf'),plot=p1,width=4,height=4)
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
  write_tsv("data/B22_summary-of-subjects.tsv")

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
  write_tsv('fig/B22_22_dx_DBS_target.tsv')

setdiff(DBS_aper$subject,DBS_aper_agg2$subject)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,x=dx_dbs_target,fill=dx_dbs_target)+
  geom_dotplot(binaxis = "y", stackdir = "center",binwidth=0.05,binposition='all') +
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_y_continuous(limits=c(2.5,4.1))+
  scale_fill_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B22_22_exponent-vs-dx.pdf',width=1.8,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,x=dx_dbs_target,color=dx_dbs_target,fill=dx_dbs_target)+
  geom_boxplot(alpha=0,width=0.6) +
  geom_point(position=position_jitter(width=0.25),colour="black",shape=21)+
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_y_continuous(limits=c(2.5,4.1))+
  scale_color_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))+
  scale_fill_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B22_22b_exponent-vs-dx.pdf',width=1.8,height=3)


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
ggsave('fig/B22_23b_log-knee-vs-dx.pdf',width=1.8,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=offset,x=dx_dbs_target,color=dx_dbs_target)+
  geom_boxplot(alpha=0,width=0.6) +
  geom_jitter(width=0.25)+
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_color_manual(values=c(color_pd_ecog,color_pd_ecog,color_et_ecog))
ggsave('fig/B22_23b_offset-vs-dx.pdf',width=1.8,height=3)

#---------

summary(lm(exponent ~ dx_dbs_target, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      3.12132    0.08206  38.039   <2e-16 ***
#   dx_dbs_target.L -0.15089    0.09049  -1.667    0.102    
# dx_dbs_target.Q  0.22693    0.17947   1.264    0.212   

oneway_test(exponent ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2, 
            distribution=approximate(nresample=999))
# Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by factor(dx_dbs_target, ordered = FALSE) (PD_STN, PD_GPi, ET_VIM)
# chi-squared = 4.3984, p-value = 0.1091

oneway_test(exponent ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2 %>% filter(dx_dbs_target %in% c("PD_STN","ET_VIM")), 
            distribution=approximate(nresample=999))

# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by factor(dx_dbs_target, ordered = FALSE) (PD_STN, ET_VIM)
# Z = 1.6113, p-value = 0.1151
# alternative hypothesis: true mu is not equal to 0

summary(lm(log_knee ~ dx_dbs_target, data=DBS_aper_agg2))
# Call:
#   lm(formula = log_knee ~ dx_dbs_target, data = DBS_aper_agg2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.0553  0.1275  0.2093  0.2632  0.6869 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.81839    0.13361   6.125 1.88e-07 ***
#   dx_dbs_target.L -0.09027    0.14735  -0.613    0.543    
# dx_dbs_target.Q  0.42395    0.29222   1.451    0.154    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6848 on 46 degrees of freedom
# Multiple R-squared:  0.05335,	Adjusted R-squared:  0.01219 
# F-statistic: 1.296 on 2 and 46 DF,  p-value: 0.2834

oneway_test(log_knee ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2, 
            distribution=approximate(nresample=999))

# Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by factor(dx_dbs_target, ordered = FALSE) (PD_STN, PD_GPi, ET_VIM)
# chi-squared = 2.5608, p-value = 0.2032

oneway_test(log_knee ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2 %>% filter(dx_dbs_target %in% c("PD_STN","ET_VIM")), 
            distribution=approximate(nresample=999))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by factor(dx_dbs_target, ordered = FALSE) (PD_STN, ET_VIM)
# Z = 0.64479, p-value = 0.5415
# alternative hypothesis: true mu is not equal to 0

summary(lm(offset ~ dx_dbs_target, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      2.17437    0.13866  15.681   <2e-16 ***
#   dx_dbs_target.L -0.05164    0.15292  -0.338    0.737    
# dx_dbs_target.Q -0.49303    0.30328  -1.626    0.111    


oneway_test(offset ~ factor(dx_dbs_target,ordered=FALSE), 
            data=DBS_aper_agg2, 
            distribution=approximate(nresample=999))
# Approximative K-Sample Fisher-Pitman Permutation Test
# 
# data:  offset by factor(dx_dbs_target, ordered = FALSE) (PD_STN, PD_GPi, ET_VIM)
# chi-squared = 2.6632, p-value = 0.2623

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
  write_tsv('fig/B22_24_aper-vs-age.tsv')

DBS_aper_agg24 %>%
  ggplot() +
  aes(x=dbs_age,y=exponent) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('fig/B22_24_exponent-vs-age.pdf',width=5,height=4)

DBS_aper_agg24 %>%
  with(cor.test(exponent,dbs_age,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  exponent and dbs_age
# S = 12957, p-value = 0.3371
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1464478 

DBS_aper_agg24 %>%
  ggplot() +
  aes(x=dbs_age,y=log_knee) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('fig/B22_24_log-knee-vs-age.pdf',width=5,height=4)

DBS_aper_agg24 %>%
  with(cor.test(log_knee,dbs_age,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  log_knee and dbs_age
# S = 13845, p-value = 0.5657
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.08793063 

DBS_aper_agg24 %>%
  ggplot() +
  aes(x=dbs_age,y=offset) +
  geom_point() +
  geom_smooth(method='lm')
ggsave('fig/B22_24_offset-vs-age.pdf',width=5,height=4)

DBS_aper_agg24 %>%
  with(cor.test(offset,dbs_age,method = 'spearman'))

# Spearman's rank correlation rho
# 
# data:  offset and dbs_age
# S = 13482, p-value = 0.4644
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1118808 


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
#   1              3.36              1.25            1.76


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
ggsave('fig/B22_28_aper_log-knee_vs_anatomy.pdf',width=5,height=4)

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
ggsave('fig/B22_28b_aper_log-knee_vs_anatomy.pdf',width=5,height=4)

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
  write_tsv('data/B22_aper_param_summary_by_MMP1.tsv')

me_mod <- lmer(log_knee ~ region + (1 | subject), 
               data=aper_ecog_roi_agg %>% within(region<-relevel(region,'6v')))
summary(me_mod)
# 
# Fixed effects:
#   Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)   1.246e+00  1.654e-02  1.553e+02  75.357   <2e-16 ***
#   regionp9-46v  2.633e-02  2.863e-02  3.273e+02   0.919   0.3585    
# region8C      1.497e-02  2.122e-02  3.253e+02   0.705   0.4812    
# region8Av    -2.121e-02  2.240e-02  3.254e+02  -0.946   0.3446    
# regionIFSp    6.114e-03  2.643e-02  3.266e+02   0.231   0.8172    
# region44      7.241e-03  2.406e-02  3.260e+02   0.301   0.7636    
# region6r     -3.820e-02  2.403e-02  3.255e+02  -1.590   0.1129    
# region55b     2.028e-02  1.998e-02  3.247e+02   1.015   0.3109    
# region4       3.619e-02  2.398e-02  3.251e+02   1.509   0.1322    
# region3a     -4.896e-02  2.633e-02  3.260e+02  -1.860   0.0638 .  
# region3b      3.219e-02  2.130e-02  3.269e+02   1.512   0.1316    
# region1       4.305e-02  1.752e-02  3.234e+02   2.457   0.0145 *  
#   region43      3.827e-03  2.338e-02  3.251e+02   0.164   0.8701    
# regionOP4    -2.613e-02  1.943e-02  3.237e+02  -1.345   0.1797    
# regionPFop   -1.013e-02  2.003e-02  3.255e+02  -0.506   0.6133    
# regionPF      3.540e-05  2.265e-02  3.293e+02   0.002   0.9988    
# regionPSL    -4.887e-02  2.880e-02  3.290e+02  -1.697   0.0907 .  
# regionSTV    -1.151e-02  2.702e-02  3.420e+02  -0.426   0.6703    
# regionA5     -5.337e-02  2.066e-02  3.277e+02  -2.584   0.0102 *  
#   regionA4     -1.125e-02  1.822e-02  3.269e+02  -0.618   0.5373    
# ---


amod <- aov(log_knee_s ~ region, data=aper_ecog_roi_agg)
amod
summary(amod)
tuk <- glht(amod, linfct = mcp(region = "Tukey"))
summary(tuk) 
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: Tukey Contrasts
# 
# 
# Fit: aov(formula = log_knee_s ~ region, data = aper_ecog_roi_agg)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
# 8C - p9-46v == 0   -0.0084210  0.0293842  -0.287   1.0000    
# 8Av - p9-46v == 0  -0.0436795  0.0301781  -1.447   0.9962    
# IFSp - p9-46v == 0 -0.0215627  0.0329040  -0.655   1.0000    
# 44 - p9-46v == 0   -0.0124261  0.0312773  -0.397   1.0000    
# 6r - p9-46v == 0   -0.0551919  0.0312773  -1.765   0.9652    
# 6v - p9-46v == 0   -0.0219015  0.0273604  -0.800   1.0000    
# 55b - p9-46v == 0  -0.0025233  0.0286141  -0.088   1.0000    
# 4 - p9-46v == 0     0.0155033  0.0312773   0.496   1.0000    
# 3a - p9-46v == 0   -0.0779205  0.0329040  -2.368   0.6708    
# 3b - p9-46v == 0    0.0053848  0.0293842   0.183   1.0000    
# 1 - p9-46v == 0     0.0234077  0.0271385   0.863   1.0000    
# 43 - p9-46v == 0   -0.0195156  0.0308666  -0.632   1.0000    
# OP4 - p9-46v == 0  -0.0455325  0.0283124  -1.608   0.9869    
# PFop - p9-46v == 0 -0.0330767  0.0286141  -1.156   0.9998    
# PF - p9-46v == 0   -0.0221630  0.0301781  -0.734   1.0000    
# PSL - p9-46v == 0  -0.0691933  0.0345100  -2.005   0.8915    
# STV - p9-46v == 0  -0.0414115  0.0329040  -1.259   0.9994    
# A5 - p9-46v == 0   -0.0710112  0.0289667  -2.451   0.6068    
# A4 - p9-46v == 0   -0.0292133  0.0275294  -1.061   0.9999    
# 8Av - 8C == 0      -0.0352585  0.0241497  -1.460   0.9957    
# IFSp - 8C == 0     -0.0131416  0.0274802  -0.478   1.0000    
# 44 - 8C == 0       -0.0040050  0.0255100  -0.157   1.0000    
# 6r - 8C == 0       -0.0467709  0.0255100  -1.833   0.9503    
# 6v - 8C == 0       -0.0134804  0.0205203  -0.657   1.0000    
# 55b - 8C == 0       0.0058977  0.0221644   0.266   1.0000    
# 4 - 8C == 0         0.0239243  0.0255100   0.938   1.0000    
# 3a - 8C == 0       -0.0694995  0.0274802  -2.529   0.5475    
# 3b - 8C == 0        0.0138059  0.0231500   0.596   1.0000    
# 1 - 8C == 0         0.0318287  0.0202236   1.574   0.9896    
# 43 - 8C == 0       -0.0110946  0.0250048  -0.444   1.0000    
# OP4 - 8C == 0      -0.0371114  0.0217735  -1.704   0.9756    
# PFop - 8C == 0     -0.0246557  0.0221644  -1.112   0.9999    
# PF - 8C == 0       -0.0137419  0.0241497  -0.569   1.0000    
# PSL - 8C == 0      -0.0607723  0.0293842  -2.068   0.8629    
# STV - 8C == 0      -0.0329904  0.0274802  -1.201   0.9997    
# A5 - 8C == 0       -0.0625902  0.0226177  -2.767   0.3676    
# A4 - 8C == 0       -0.0207923  0.0207452  -1.002   1.0000    
# IFSp - 8Av == 0     0.0221169  0.0283275   0.781   1.0000    
# 44 - 8Av == 0       0.0312535  0.0264206   1.183   0.9997    
# 6r - 8Av == 0      -0.0115124  0.0264206  -0.436   1.0000    
# 6v - 8Av == 0       0.0217780  0.0216418   1.006   1.0000    
# 55b - 8Av == 0      0.0411562  0.0232067   1.773   0.9637    
# 4 - 8Av == 0        0.0591828  0.0264206   2.240   0.7616    
# 3a - 8Av == 0      -0.0342410  0.0283275  -1.209   0.9996    
# 3b - 8Av == 0       0.0490644  0.0241497   2.032   0.8796    
# 1 - 8Av == 0        0.0670872  0.0213607   3.141   0.1613    
# 43 - 8Av == 0       0.0241639  0.0259332   0.932   1.0000    
# OP4 - 8Av == 0     -0.0018529  0.0228336  -0.081   1.0000    
# PFop - 8Av == 0     0.0106028  0.0232067   0.457   1.0000    
# PF - 8Av == 0       0.0215166  0.0251097   0.857   1.0000    
# PSL - 8Av == 0     -0.0255138  0.0301781  -0.845   1.0000    
# STV - 8Av == 0      0.0022680  0.0283275   0.080   1.0000    
# A5 - 8Av == 0      -0.0273317  0.0236400  -1.156   0.9998    
# A4 - 8Av == 0       0.0144662  0.0218552   0.662   1.0000    
# 44 - IFSp == 0      0.0091366  0.0294958   0.310   1.0000    
# 6r - IFSp == 0     -0.0336293  0.0294958  -1.140   0.9998    
# 6v - IFSp == 0     -0.0003388  0.0253046  -0.013   1.0000    
# 55b - IFSp == 0     0.0190393  0.0266552   0.714   1.0000    
# 4 - IFSp == 0       0.0370659  0.0294958   1.257   0.9994    
# 3a - IFSp == 0     -0.0563578  0.0312154  -1.805   0.9565    
# 3b - IFSp == 0      0.0269475  0.0274802   0.981   1.0000    
# 1 - IFSp == 0       0.0449703  0.0250646   1.794   0.9584    
# 43 - IFSp == 0      0.0020470  0.0290600   0.070   1.0000    
# OP4 - IFSp == 0    -0.0239698  0.0263311  -0.910   1.0000    
# PFop - IFSp == 0   -0.0115141  0.0266552  -0.432   1.0000    
# PF - IFSp == 0     -0.0006003  0.0283275  -0.021   1.0000    
# PSL - IFSp == 0    -0.0476307  0.0329040  -1.448   0.9962    
# STV - IFSp == 0    -0.0198488  0.0312154  -0.636   1.0000    
# A5 - IFSp == 0     -0.0494485  0.0270334  -1.829   0.9515    
# A4 - IFSp == 0     -0.0076507  0.0254873  -0.300   1.0000    
# 6r - 44 == 0       -0.0427658  0.0276695  -1.546   0.9917    
# 6v - 44 == 0       -0.0094754  0.0231500  -0.409   1.0000    
# 55b - 44 == 0       0.0099027  0.0246191   0.402   1.0000    
# 4 - 44 == 0         0.0279293  0.0276695   1.009   1.0000    
# 3a - 44 == 0       -0.0654944  0.0294958  -2.220   0.7743    
# 3b - 44 == 0        0.0178109  0.0255100   0.698   1.0000    
# 1 - 44 == 0         0.0358338  0.0228874   1.566   0.9903    
# 43 - 44 == 0       -0.0070896  0.0272045  -0.261   1.0000    
# OP4 - 44 == 0      -0.0331064  0.0242678  -1.364   0.9982    
# PFop - 44 == 0     -0.0206507  0.0246191  -0.839   1.0000    
# PF - 44 == 0       -0.0097369  0.0264206  -0.369   1.0000    
# PSL - 44 == 0      -0.0567672  0.0312773  -1.815   0.9544    
# STV - 44 == 0      -0.0289854  0.0294958  -0.983   1.0000    
# A5 - 44 == 0       -0.0585851  0.0250280  -2.341   0.6899    
# A4 - 44 == 0       -0.0167873  0.0233496  -0.719   1.0000    
# 6v - 6r == 0        0.0332904  0.0231500   1.438   0.9965    
# 55b - 6r == 0       0.0526686  0.0246191   2.139   0.8253    
# 4 - 6r == 0         0.0706952  0.0276695   2.555   0.5267    
# 3a - 6r == 0       -0.0227286  0.0294958  -0.771   1.0000    
# 3b - 6r == 0        0.0605767  0.0255100   2.375   0.6639    
# 1 - 6r == 0         0.0785996  0.0228874   3.434   0.0708 .  
# 43 - 6r == 0        0.0356763  0.0272045   1.311   0.9990    
# OP4 - 6r == 0       0.0096594  0.0242678   0.398   1.0000    
# PFop - 6r == 0      0.0221152  0.0246191   0.898   1.0000    
# PF - 6r == 0        0.0330289  0.0264206   1.250   0.9994    
# PSL - 6r == 0      -0.0140014  0.0312773  -0.448   1.0000    
# STV - 6r == 0       0.0137804  0.0294958   0.467   1.0000    
# A5 - 6r == 0       -0.0158193  0.0250280  -0.632   1.0000    
# A4 - 6r == 0        0.0259786  0.0233496   1.113   0.9999    
# 55b - 6v == 0       0.0193781  0.0194016   0.999   1.0000    
# 4 - 6v == 0         0.0374047  0.0231500   1.616   0.9862    
# 3a - 6v == 0       -0.0560190  0.0253046  -2.214   0.7783    
# 3b - 6v == 0        0.0272863  0.0205203   1.330   0.9987    
# 1 - 6v == 0         0.0453092  0.0171509   2.642   0.4609    
# 43 - 6v == 0        0.0023858  0.0225921   0.106   1.0000    
# OP4 - 6v == 0      -0.0236310  0.0189537  -1.247   0.9995    
# PFop - 6v == 0     -0.0111753  0.0194016  -0.576   1.0000    
# PF - 6v == 0       -0.0002615  0.0216418  -0.012   1.0000    
# PSL - 6v == 0      -0.0472918  0.0273604  -1.728   0.9716    
# STV - 6v == 0      -0.0195100  0.0253046  -0.771   1.0000    
# A5 - 6v == 0       -0.0491097  0.0199179  -2.466   0.5970    
# A4 - 6v == 0       -0.0073118  0.0177629  -0.412   1.0000    
# 4 - 55b == 0        0.0180266  0.0246191   0.732   1.0000    
# 3a - 55b == 0      -0.0753971  0.0266552  -2.829   0.3273    
# 3b - 55b == 0       0.0079082  0.0221644   0.357   1.0000    
# 1 - 55b == 0        0.0259310  0.0190875   1.359   0.9983    
# 43 - 55b == 0      -0.0169923  0.0240953  -0.705   1.0000    
# OP4 - 55b == 0     -0.0430091  0.0207226  -2.075   0.8588    
# PFop - 55b == 0    -0.0305534  0.0211329  -1.446   0.9963    
# PF - 55b == 0      -0.0196396  0.0232067  -0.846   1.0000    
# PSL - 55b == 0     -0.0666700  0.0286141  -2.330   0.6984    
# STV - 55b == 0     -0.0388881  0.0266552  -1.459   0.9959    
# A5 - 55b == 0      -0.0684878  0.0216079  -3.170   0.1491    
# A4 - 55b == 0      -0.0266900  0.0196393  -1.359   0.9983    
# 3a - 4 == 0        -0.0934238  0.0294958  -3.167   0.1506    
# 3b - 4 == 0        -0.0101184  0.0255100  -0.397   1.0000    
# 1 - 4 == 0          0.0079044  0.0228874   0.345   1.0000    
# 43 - 4 == 0        -0.0350189  0.0272045  -1.287   0.9992    
# OP4 - 4 == 0       -0.0610357  0.0242678  -2.515   0.5574    
# PFop - 4 == 0      -0.0485800  0.0246191  -1.973   0.9043    
# PF - 4 == 0        -0.0376662  0.0264206  -1.426   0.9968    
# PSL - 4 == 0       -0.0846966  0.0312773  -2.708   0.4091    
# STV - 4 == 0       -0.0569148  0.0294958  -1.930   0.9206    
# A5 - 4 == 0        -0.0865145  0.0250280  -3.457   0.0661 .  
# A4 - 4 == 0        -0.0447166  0.0233496  -1.915   0.9261    
# 3b - 3a == 0        0.0833053  0.0274802   3.031   0.2093    
# 1 - 3a == 0         0.1013282  0.0250646   4.043    <0.01 ** 
#   43 - 3a == 0        0.0584048  0.0290600   2.010   0.8894    
# OP4 - 3a == 0       0.0323880  0.0263311   1.230   0.9995    
# PFop - 3a == 0      0.0448437  0.0266552   1.682   0.9787    
# PF - 3a == 0        0.0557575  0.0283275   1.968   0.9068    
# PSL - 3a == 0       0.0087272  0.0329040   0.265   1.0000    
# STV - 3a == 0       0.0365090  0.0312154   1.170   0.9998    
# A5 - 3a == 0        0.0069093  0.0270334   0.256   1.0000    
# A4 - 3a == 0        0.0487072  0.0254873   1.911   0.9274    
# 1 - 3b == 0         0.0180228  0.0202236   0.891   1.0000    
# 43 - 3b == 0       -0.0249005  0.0250048  -0.996   1.0000    
# OP4 - 3b == 0      -0.0509173  0.0217735  -2.339   0.6925    
# PFop - 3b == 0     -0.0384616  0.0221644  -1.735   0.9708    
# PF - 3b == 0       -0.0275478  0.0241497  -1.141   0.9998    
# PSL - 3b == 0      -0.0745782  0.0293842  -2.538   0.5412    
# STV - 3b == 0      -0.0467963  0.0274802  -1.703   0.9757    
# A5 - 3b == 0       -0.0763960  0.0226177  -3.378   0.0851 .  
# A4 - 3b == 0       -0.0345982  0.0207452  -1.668   0.9805    
# 43 - 1 == 0        -0.0429233  0.0223229  -1.923   0.9230    
# OP4 - 1 == 0       -0.0689402  0.0186321  -3.700   0.0304 *  
#   PFop - 1 == 0      -0.0564844  0.0190875  -2.959   0.2479    
# PF - 1 == 0        -0.0455706  0.0213607  -2.133   0.8271    
# PSL - 1 == 0       -0.0926010  0.0271385  -3.412   0.0762 .  
# STV - 1 == 0       -0.0648192  0.0250646  -2.586   0.5032    
# A5 - 1 == 0        -0.0944189  0.0196121  -4.814    <0.01 ***
#   A4 - 1 == 0        -0.0526210  0.0174193  -3.021   0.2152    
# OP4 - 43 == 0      -0.0260168  0.0237362  -1.096   0.9999    
# PFop - 43 == 0     -0.0135611  0.0240953  -0.563   1.0000    
# PF - 43 == 0       -0.0026473  0.0259332  -0.102   1.0000    
# PSL - 43 == 0      -0.0496777  0.0308666  -1.609   0.9868    
# STV - 43 == 0      -0.0218958  0.0290600  -0.753   1.0000    
# A5 - 43 == 0       -0.0514956  0.0245129  -2.101   0.8456    
# A4 - 43 == 0       -0.0096977  0.0227965  -0.425   1.0000    
# PFop - OP4 == 0     0.0124557  0.0207226   0.601   1.0000    
# PF - OP4 == 0       0.0233695  0.0228336   1.023   1.0000    
# PSL - OP4 == 0     -0.0236608  0.0283124  -0.836   1.0000    
# STV - OP4 == 0      0.0041210  0.0263311   0.157   1.0000    
# A5 - OP4 == 0      -0.0254787  0.0212067  -1.201   0.9997    
# A4 - OP4 == 0       0.0163191  0.0191970   0.850   1.0000    
# PF - PFop == 0      0.0109138  0.0232067   0.470   1.0000    
# PSL - PFop == 0    -0.0361166  0.0286141  -1.262   0.9994    
# STV - PFop == 0    -0.0083347  0.0266552  -0.313   1.0000    
# A5 - PFop == 0     -0.0379345  0.0216079  -1.756   0.9671    
# A4 - PFop == 0      0.0038634  0.0196393   0.197   1.0000    
# PSL - PF == 0      -0.0470304  0.0301781  -1.558   0.9909    
# STV - PF == 0      -0.0192485  0.0283275  -0.679   1.0000    
# A5 - PF == 0       -0.0488482  0.0236400  -2.066   0.8627    
# A4 - PF == 0       -0.0070504  0.0218552  -0.323   1.0000    
# STV - PSL == 0      0.0277818  0.0329040   0.844   1.0000    
# A5 - PSL == 0      -0.0018179  0.0289667  -0.063   1.0000    
# A4 - PSL == 0       0.0399800  0.0275294   1.452   0.9960    
# A5 - STV == 0      -0.0295997  0.0270334  -1.095   0.9999    
# A4 - STV == 0       0.0121982  0.0254873   0.479   1.0000    
# A4 - A5 == 0        0.0417979  0.0201495   2.074   0.8588    
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
ggsave('fig/B22_28c_aper_exponent_vs_anatomy.pdf',width=5,height=4)

amod <- aov(exponent_s ~ region, data=aper_ecog_roi_agg)
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
# Fit: aov(formula = exponent_s ~ region, data = aper_ecog_roi_agg)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)   
# 8C - p9-46v == 0   -5.872e-02  1.043e-01  -0.563   1.0000   
# 8Av - p9-46v == 0  -1.976e-01  1.071e-01  -1.845   0.9470   
# IFSp - p9-46v == 0 -1.701e-01  1.168e-01  -1.457   0.9958   
# 44 - p9-46v == 0   -2.564e-01  1.110e-01  -2.311   0.7127   
# 6r - p9-46v == 0   -3.147e-01  1.110e-01  -2.835   0.3232   
# 6v - p9-46v == 0   -3.107e-03  9.708e-02  -0.032   1.0000   
# 55b - p9-46v == 0  -1.237e-02  1.015e-01  -0.122   1.0000   
# 4 - p9-46v == 0    -8.213e-02  1.110e-01  -0.740   1.0000   
# 3a - p9-46v == 0   -2.540e-01  1.168e-01  -2.175   0.8037   
# 3b - p9-46v == 0   -4.548e-02  1.043e-01  -0.436   1.0000   
# 1 - p9-46v == 0    -7.269e-02  9.630e-02  -0.755   1.0000   
# 43 - p9-46v == 0   -2.529e-01  1.095e-01  -2.309   0.7141   
# OP4 - p9-46v == 0  -2.940e-01  1.005e-01  -2.927   0.2661   
# PFop - p9-46v == 0 -2.158e-01  1.015e-01  -2.125   0.8315   
# PF - p9-46v == 0   -2.703e-01  1.071e-01  -2.525   0.5482   
# PSL - p9-46v == 0  -3.462e-01  1.225e-01  -2.827   0.3277   
# STV - p9-46v == 0  -2.444e-01  1.168e-01  -2.094   0.8488   
# A5 - p9-46v == 0   -2.529e-01  1.028e-01  -2.461   0.6009   
# A4 - p9-46v == 0   -2.628e-01  9.768e-02  -2.690   0.4237   
# 8Av - 8C == 0      -1.389e-01  8.569e-02  -1.620   0.9856   
# IFSp - 8C == 0     -1.114e-01  9.751e-02  -1.142   0.9998   
# 44 - 8C == 0       -1.977e-01  9.052e-02  -2.184   0.7970   
# 6r - 8C == 0       -2.560e-01  9.052e-02  -2.828   0.3280   
# 6v - 8C == 0        5.561e-02  7.281e-02   0.764   1.0000   
# 55b - 8C == 0       4.635e-02  7.865e-02   0.589   1.0000   
# 4 - 8C == 0        -2.341e-02  9.052e-02  -0.259   1.0000   
# 3a - 8C == 0       -1.952e-01  9.751e-02  -2.002   0.8915   
# 3b - 8C == 0        1.324e-02  8.214e-02   0.161   1.0000   
# 1 - 8C == 0        -1.397e-02  7.176e-02  -0.195   1.0000   
# 43 - 8C == 0       -1.942e-01  8.873e-02  -2.189   0.7931   
# OP4 - 8C == 0      -2.353e-01  7.726e-02  -3.046   0.2038   
# PFop - 8C == 0     -1.571e-01  7.865e-02  -1.997   0.8946   
# PF - 8C == 0       -2.116e-01  8.569e-02  -2.470   0.5943   
# PSL - 8C == 0      -2.875e-01  1.043e-01  -2.757   0.3747   
# STV - 8C == 0      -1.857e-01  9.751e-02  -1.905   0.9296   
# A5 - 8C == 0       -1.942e-01  8.026e-02  -2.420   0.6308   
# A4 - 8C == 0       -2.041e-01  7.361e-02  -2.772   0.3660   
# IFSp - 8Av == 0     2.746e-02  1.005e-01   0.273   1.0000   
# 44 - 8Av == 0      -5.887e-02  9.375e-02  -0.628   1.0000   
# 6r - 8Av == 0      -1.171e-01  9.375e-02  -1.249   0.9994   
# 6v - 8Av == 0       1.945e-01  7.679e-02   2.532   0.5453   
# 55b - 8Av == 0      1.852e-01  8.235e-02   2.249   0.7555   
# 4 - 8Av == 0        1.154e-01  9.375e-02   1.231   0.9995   
# 3a - 8Av == 0      -5.639e-02  1.005e-01  -0.561   1.0000   
# 3b - 8Av == 0       1.521e-01  8.569e-02   1.775   0.9633   
# 1 - 8Av == 0        1.249e-01  7.580e-02   1.648   0.9830   
# 43 - 8Av == 0      -5.537e-02  9.202e-02  -0.602   1.0000   
# OP4 - 8Av == 0     -9.646e-02  8.102e-02  -1.191   0.9997   
# PFop - 8Av == 0    -1.820e-02  8.235e-02  -0.221   1.0000   
# PF - 8Av == 0      -7.276e-02  8.910e-02  -0.817   1.0000   
# PSL - 8Av == 0     -1.486e-01  1.071e-01  -1.388   0.9977   
# STV - 8Av == 0     -4.687e-02  1.005e-01  -0.466   1.0000   
# A5 - 8Av == 0      -5.533e-02  8.388e-02  -0.660   1.0000   
# A4 - 8Av == 0      -6.520e-02  7.755e-02  -0.841   1.0000   
# 44 - IFSp == 0     -8.633e-02  1.047e-01  -0.825   1.0000   
# 6r - IFSp == 0     -1.446e-01  1.047e-01  -1.381   0.9979   
# 6v - IFSp == 0      1.670e-01  8.979e-02   1.860   0.9427   
# 55b - IFSp == 0     1.577e-01  9.458e-02   1.668   0.9805   
# 4 - IFSp == 0       8.799e-02  1.047e-01   0.841   1.0000   
# 3a - IFSp == 0     -8.385e-02  1.108e-01  -0.757   1.0000   
# 3b - IFSp == 0      1.246e-01  9.751e-02   1.278   0.9992   
# 1 - IFSp == 0       9.743e-02  8.894e-02   1.095   0.9999   
# 43 - IFSp == 0     -8.283e-02  1.031e-01  -0.803   1.0000   
# OP4 - IFSp == 0    -1.239e-01  9.343e-02  -1.326   0.9987   
# PFop - IFSp == 0   -4.566e-02  9.458e-02  -0.483   1.0000   
# PF - IFSp == 0     -1.002e-01  1.005e-01  -0.997   1.0000   
# PSL - IFSp == 0    -1.761e-01  1.168e-01  -1.508   0.9937   
# STV - IFSp == 0    -7.433e-02  1.108e-01  -0.671   1.0000   
# A5 - IFSp == 0     -8.279e-02  9.592e-02  -0.863   1.0000   
# A4 - IFSp == 0     -9.266e-02  9.044e-02  -1.025   1.0000   
# 6r - 44 == 0       -5.823e-02  9.818e-02  -0.593   1.0000   
# 6v - 44 == 0        2.533e-01  8.214e-02   3.084   0.1865   
# 55b - 44 == 0       2.441e-01  8.736e-02   2.794   0.3504   
# 4 - 44 == 0         1.743e-01  9.818e-02   1.775   0.9631   
# 3a - 44 == 0        2.477e-03  1.047e-01   0.024   1.0000   
# 3b - 44 == 0        2.110e-01  9.052e-02   2.331   0.6988   
# 1 - 44 == 0         1.838e-01  8.121e-02   2.263   0.7469   
# 43 - 44 == 0        3.498e-03  9.653e-02   0.036   1.0000   
# OP4 - 44 == 0      -3.759e-02  8.611e-02  -0.437   1.0000   
# PFop - 44 == 0      4.067e-02  8.736e-02   0.466   1.0000   
# PF - 44 == 0       -1.389e-02  9.375e-02  -0.148   1.0000   
# PSL - 44 == 0      -8.975e-02  1.110e-01  -0.809   1.0000   
# STV - 44 == 0       1.200e-02  1.047e-01   0.115   1.0000   
# A5 - 44 == 0        3.539e-03  8.881e-02   0.040   1.0000   
# A4 - 44 == 0       -6.335e-03  8.285e-02  -0.076   1.0000   
# 6v - 6r == 0        3.116e-01  8.214e-02   3.793   0.0226 * 
#   55b - 6r == 0       3.023e-01  8.736e-02   3.460   0.0670 . 
# 4 - 6r == 0         2.325e-01  9.818e-02   2.369   0.6696   
# 3a - 6r == 0        6.071e-02  1.047e-01   0.580   1.0000   
# 3b - 6r == 0        2.692e-01  9.052e-02   2.974   0.2401   
# 1 - 6r == 0         2.420e-01  8.121e-02   2.980   0.2368   
# 43 - 6r == 0        6.173e-02  9.653e-02   0.639   1.0000   
# OP4 - 6r == 0       2.064e-02  8.611e-02   0.240   1.0000   
# PFop - 6r == 0      9.890e-02  8.736e-02   1.132   0.9999   
# PF - 6r == 0        4.434e-02  9.375e-02   0.473   1.0000   
# PSL - 6r == 0      -3.153e-02  1.110e-01  -0.284   1.0000   
# STV - 6r == 0       7.023e-02  1.047e-01   0.671   1.0000   
# A5 - 6r == 0        6.177e-02  8.881e-02   0.696   1.0000   
# A4 - 6r == 0        5.189e-02  8.285e-02   0.626   1.0000   
# 55b - 6v == 0      -9.268e-03  6.884e-02  -0.135   1.0000   
# 4 - 6v == 0        -7.902e-02  8.214e-02  -0.962   1.0000   
# 3a - 6v == 0       -2.509e-01  8.979e-02  -2.794   0.3519   
# 3b - 6v == 0       -4.237e-02  7.281e-02  -0.582   1.0000   
# 1 - 6v == 0        -6.958e-02  6.086e-02  -1.143   0.9998   
# 43 - 6v == 0       -2.498e-01  8.016e-02  -3.117   0.1711   
# OP4 - 6v == 0      -2.909e-01  6.725e-02  -4.326    <0.01 **
#   PFop - 6v == 0     -2.127e-01  6.884e-02  -3.089   0.1833   
# PF - 6v == 0       -2.672e-01  7.679e-02  -3.480   0.0629 . 
# PSL - 6v == 0      -3.431e-01  9.708e-02  -3.534   0.0538 . 
# STV - 6v == 0      -2.413e-01  8.979e-02  -2.688   0.4259   
# A5 - 6v == 0       -2.498e-01  7.068e-02  -3.534   0.0527 . 
# A4 - 6v == 0       -2.597e-01  6.303e-02  -4.120    <0.01 **
#   4 - 55b == 0       -6.976e-02  8.736e-02  -0.799   1.0000   
# 3a - 55b == 0      -2.416e-01  9.458e-02  -2.554   0.5290   
# 3b - 55b == 0      -3.311e-02  7.865e-02  -0.421   1.0000   
# 1 - 55b == 0       -6.031e-02  6.773e-02  -0.891   1.0000   
# 43 - 55b == 0      -2.406e-01  8.550e-02  -2.814   0.3362   
# OP4 - 55b == 0     -2.817e-01  7.353e-02  -3.831   0.0207 * 
#   PFop - 55b == 0    -2.034e-01  7.499e-02  -2.713   0.4062   
# PF - 55b == 0      -2.580e-01  8.235e-02  -3.133   0.1649   
# PSL - 55b == 0     -3.338e-01  1.015e-01  -3.288   0.1076   
# STV - 55b == 0     -2.321e-01  9.458e-02  -2.454   0.6056   
# A5 - 55b == 0      -2.405e-01  7.667e-02  -3.137   0.1642   
# A4 - 55b == 0      -2.504e-01  6.969e-02  -3.593   0.0439 * 
#   3a - 4 == 0        -1.718e-01  1.047e-01  -1.642   0.9837   
# 3b - 4 == 0         3.665e-02  9.052e-02   0.405   1.0000   
# 1 - 4 == 0          9.441e-03  8.121e-02   0.116   1.0000   
# 43 - 4 == 0        -1.708e-01  9.653e-02  -1.770   0.9641   
# OP4 - 4 == 0       -2.119e-01  8.611e-02  -2.461   0.5998   
# PFop - 4 == 0      -1.336e-01  8.736e-02  -1.530   0.9926   
# PF - 4 == 0        -1.882e-01  9.375e-02  -2.008   0.8910   
# PSL - 4 == 0       -2.641e-01  1.110e-01  -2.379   0.6621   
# STV - 4 == 0       -1.623e-01  1.047e-01  -1.551   0.9914   
# A5 - 4 == 0        -1.708e-01  8.881e-02  -1.923   0.9239   
# A4 - 4 == 0        -1.807e-01  8.285e-02  -2.180   0.7981   
# 3b - 3a == 0        2.085e-01  9.751e-02   2.138   0.8249   
# 1 - 3a == 0         1.813e-01  8.894e-02   2.038   0.8768   
# 43 - 3a == 0        1.020e-03  1.031e-01   0.010   1.0000   
# OP4 - 3a == 0      -4.007e-02  9.343e-02  -0.429   1.0000   
# PFop - 3a == 0      3.819e-02  9.458e-02   0.404   1.0000   
# PF - 3a == 0       -1.637e-02  1.005e-01  -0.163   1.0000   
# PSL - 3a == 0      -9.223e-02  1.168e-01  -0.790   1.0000   
# STV - 3a == 0       9.519e-03  1.108e-01   0.086   1.0000   
# A5 - 3a == 0        1.062e-03  9.592e-02   0.011   1.0000   
# A4 - 3a == 0       -8.813e-03  9.044e-02  -0.097   1.0000   
# 1 - 3b == 0        -2.721e-02  7.176e-02  -0.379   1.0000   
# 43 - 3b == 0       -2.075e-01  8.873e-02  -2.338   0.6917   
# OP4 - 3b == 0      -2.486e-01  7.726e-02  -3.217   0.1321   
# PFop - 3b == 0     -1.703e-01  7.865e-02  -2.165   0.8082   
# PF - 3b == 0       -2.249e-01  8.569e-02  -2.624   0.4738   
# PSL - 3b == 0      -3.007e-01  1.043e-01  -2.884   0.2926   
# STV - 3b == 0      -1.990e-01  9.751e-02  -2.041   0.8758   
# A5 - 3b == 0       -2.074e-01  8.026e-02  -2.585   0.5041   
# A4 - 3b == 0       -2.173e-01  7.361e-02  -2.952   0.2524   
# 43 - 1 == 0        -1.803e-01  7.921e-02  -2.276   0.7368   
# OP4 - 1 == 0       -2.213e-01  6.611e-02  -3.348   0.0916 . 
# PFop - 1 == 0      -1.431e-01  6.773e-02  -2.113   0.8381   
# PF - 1 == 0        -1.976e-01  7.580e-02  -2.608   0.4867   
# PSL - 1 == 0       -2.735e-01  9.630e-02  -2.840   0.3187   
# STV - 1 == 0       -1.718e-01  8.894e-02  -1.931   0.9206   
# A5 - 1 == 0        -1.802e-01  6.959e-02  -2.590   0.5000   
# A4 - 1 == 0        -1.901e-01  6.181e-02  -3.075   0.1887   
# OP4 - 43 == 0      -4.109e-02  8.422e-02  -0.488   1.0000   
# PFop - 43 == 0      3.717e-02  8.550e-02   0.435   1.0000   
# PF - 43 == 0       -1.739e-02  9.202e-02  -0.189   1.0000   
# PSL - 43 == 0      -9.325e-02  1.095e-01  -0.851   1.0000   
# STV - 43 == 0       8.499e-03  1.031e-01   0.082   1.0000   
# A5 - 43 == 0        4.183e-05  8.698e-02   0.000   1.0000   
# A4 - 43 == 0       -9.833e-03  8.089e-02  -0.122   1.0000   
# PFop - OP4 == 0     7.826e-02  7.353e-02   1.064   0.9999   
# PF - OP4 == 0       2.370e-02  8.102e-02   0.293   1.0000   
# PSL - OP4 == 0     -5.216e-02  1.005e-01  -0.519   1.0000   
# STV - OP4 == 0      4.959e-02  9.343e-02   0.531   1.0000   
# A5 - OP4 == 0       4.113e-02  7.525e-02   0.547   1.0000   
# A4 - OP4 == 0       3.125e-02  6.812e-02   0.459   1.0000   
# PF - PFop == 0     -5.456e-02  8.235e-02  -0.663   1.0000   
# PSL - PFop == 0    -1.304e-01  1.015e-01  -1.285   0.9992   
# STV - PFop == 0    -2.867e-02  9.458e-02  -0.303   1.0000   
# A5 - PFop == 0     -3.713e-02  7.667e-02  -0.484   1.0000   
# A4 - PFop == 0     -4.700e-02  6.969e-02  -0.675   1.0000   
# PSL - PF == 0      -7.586e-02  1.071e-01  -0.708   1.0000   
# STV - PF == 0       2.589e-02  1.005e-01   0.258   1.0000   
# A5 - PF == 0        1.743e-02  8.388e-02   0.208   1.0000   
# A4 - PF == 0        7.555e-03  7.755e-02   0.097   1.0000   
# STV - PSL == 0      1.018e-01  1.168e-01   0.871   1.0000   
# A5 - PSL == 0       9.329e-02  1.028e-01   0.908   1.0000   
# A4 - PSL == 0       8.342e-02  9.768e-02   0.854   1.0000   
# A5 - STV == 0      -8.457e-03  9.592e-02  -0.088   1.0000   
# A4 - STV == 0      -1.833e-02  9.044e-02  -0.203   1.0000   
# A4 - A5 == 0       -9.875e-03  7.150e-02  -0.138   1.0000   
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
ggsave('fig/B22_28c_aper_offset_vs_anatomy.pdf',width=5,height=4)

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
# 8C - p9-46v == 0   -0.134053   0.115592  -1.160   0.9998    
# 8Av - p9-46v == 0  -0.071489   0.118715  -0.602   1.0000    
# IFSp - p9-46v == 0 -0.166870   0.129438  -1.289   0.9992    
# 44 - p9-46v == 0   -0.391743   0.123039  -3.184   0.1443    
# 6r - p9-46v == 0   -0.438185   0.123039  -3.561   0.0478 *  
#   6v - p9-46v == 0   -0.005639   0.107631  -0.052   1.0000    
# 55b - p9-46v == 0  -0.028062   0.112563  -0.249   1.0000    
# 4 - p9-46v == 0    -0.201054   0.123039  -1.634   0.9844    
# 3a - p9-46v == 0   -0.361884   0.129438  -2.796   0.3476    
# 3b - p9-46v == 0   -0.124148   0.115592  -1.074   0.9999    
# 1 - p9-46v == 0    -0.216087   0.106758  -2.024   0.8836    
# 43 - p9-46v == 0   -0.260566   0.121424  -2.146   0.8195    
# OP4 - p9-46v == 0  -0.352391   0.111376  -3.164   0.1517    
# PFop - p9-46v == 0 -0.226109   0.112563  -2.009   0.8896    
# PF - p9-46v == 0   -0.212773   0.118715  -1.792   0.9593    
# PSL - p9-46v == 0  -0.171254   0.135756  -1.261   0.9994    
# STV - p9-46v == 0  -0.173908   0.129438  -1.344   0.9985    
# A5 - p9-46v == 0   -0.141453   0.113950  -1.241   0.9995    
# A4 - p9-46v == 0   -0.205516   0.108296  -1.898   0.9318    
# 8Av - 8C == 0       0.062564   0.095001   0.659   1.0000    
# IFSp - 8C == 0     -0.032817   0.108102  -0.304   1.0000    
# 44 - 8C == 0       -0.257690   0.100352  -2.568   0.5159    
# 6r - 8C == 0       -0.304132   0.100352  -3.031   0.2122    
# 6v - 8C == 0        0.128414   0.080723   1.591   0.9884    
# 55b - 8C == 0       0.105990   0.087191   1.216   0.9996    
# 4 - 8C == 0        -0.067001   0.100352  -0.668   1.0000    
# 3a - 8C == 0       -0.227831   0.108102  -2.108   0.8412    
# 3b - 8C == 0        0.009905   0.091068   0.109   1.0000    
# 1 - 8C == 0        -0.082034   0.079556  -1.031   1.0000    
# 43 - 8C == 0       -0.126513   0.098364  -1.286   0.9992    
# OP4 - 8C == 0      -0.218338   0.085653  -2.549   0.5315    
# PFop - 8C == 0     -0.092056   0.087191  -1.056   1.0000    
# PF - 8C == 0       -0.078720   0.095001  -0.829   1.0000    
# PSL - 8C == 0      -0.037201   0.115592  -0.322   1.0000    
# STV - 8C == 0      -0.039855   0.108102  -0.369   1.0000    
# A5 - 8C == 0       -0.007400   0.088974  -0.083   1.0000    
# A4 - 8C == 0       -0.071463   0.081608  -0.876   1.0000    
# IFSp - 8Av == 0    -0.095381   0.111435  -0.856   1.0000    
# 44 - 8Av == 0      -0.320254   0.103934  -3.081   0.1864    
# 6r - 8Av == 0      -0.366696   0.103934  -3.528   0.0548 .  
# 6v - 8Av == 0       0.065850   0.085135   0.773   1.0000    
# 55b - 8Av == 0      0.043426   0.091291   0.476   1.0000    
# 4 - 8Av == 0       -0.129566   0.103934  -1.247   0.9995    
# 3a - 8Av == 0      -0.290396   0.111435  -2.606   0.4869    
# 3b - 8Av == 0      -0.052660   0.095001  -0.554   1.0000    
# 1 - 8Av == 0       -0.144598   0.084029  -1.721   0.9730    
# 43 - 8Av == 0      -0.189077   0.102016  -1.853   0.9453    
# OP4 - 8Av == 0     -0.280902   0.089823  -3.127   0.1678    
# PFop - 8Av == 0    -0.154620   0.091291  -1.694   0.9774    
# PF - 8Av == 0      -0.141285   0.098777  -1.430   0.9967    
# PSL - 8Av == 0     -0.099765   0.118715  -0.840   1.0000    
# STV - 8Av == 0     -0.102420   0.111435  -0.919   1.0000    
# A5 - 8Av == 0      -0.069965   0.092995  -0.752   1.0000    
# A4 - 8Av == 0      -0.134028   0.085974  -1.559   0.9908    
# 44 - IFSp == 0     -0.224873   0.116031  -1.938   0.9174    
# 6r - IFSp == 0     -0.271315   0.116031  -2.338   0.6917    
# 6v - IFSp == 0      0.161231   0.099544   1.620   0.9860    
# 55b - IFSp == 0     0.138808   0.104857   1.324   0.9988    
# 4 - IFSp == 0      -0.034184   0.116031  -0.295   1.0000    
# 3a - IFSp == 0     -0.195014   0.122796  -1.588   0.9886    
# 3b - IFSp == 0      0.042722   0.108102   0.395   1.0000    
# 1 - IFSp == 0      -0.049217   0.098600  -0.499   1.0000    
# 43 - IFSp == 0     -0.093696   0.114317  -0.820   1.0000    
# OP4 - IFSp == 0    -0.185521   0.103582  -1.791   0.9600    
# PFop - IFSp == 0   -0.059239   0.104857  -0.565   1.0000    
# PF - IFSp == 0     -0.045903   0.111435  -0.412   1.0000    
# PSL - IFSp == 0    -0.004384   0.129438  -0.034   1.0000    
# STV - IFSp == 0    -0.007038   0.122796  -0.057   1.0000    
# A5 - IFSp == 0      0.025417   0.106344   0.239   1.0000    
# A4 - IFSp == 0     -0.038646   0.100262  -0.385   1.0000    
# 6r - 44 == 0       -0.046442   0.108847  -0.427   1.0000    
# 6v - 44 == 0        0.386104   0.091068   4.240    <0.01 ** 
#   55b - 44 == 0       0.363680   0.096847   3.755   0.0251 *  
#   4 - 44 == 0         0.190688   0.108847   1.752   0.9676    
# 3a - 44 == 0        0.029859   0.116031   0.257   1.0000    
# 3b - 44 == 0        0.267594   0.100352   2.667   0.4428    
# 1 - 44 == 0         0.175656   0.090035   1.951   0.9129    
# 43 - 44 == 0        0.131177   0.107017   1.226   0.9996    
# OP4 - 44 == 0       0.039352   0.095465   0.412   1.0000    
# PFop - 44 == 0      0.165634   0.096847   1.710   0.9746    
# PF - 44 == 0        0.178969   0.103934   1.722   0.9728    
# PSL - 44 == 0       0.220489   0.123039   1.792   0.9594    
# STV - 44 == 0       0.217835   0.116031   1.877   0.9379    
# A5 - 44 == 0        0.250289   0.098456   2.542   0.5356    
# A4 - 44 == 0        0.186227   0.091853   2.027   0.8813    
# 6v - 6r == 0        0.432546   0.091068   4.750    <0.01 ***
#   55b - 6r == 0       0.410122   0.096847   4.235    <0.01 ** 
#   4 - 6r == 0         0.237130   0.108847   2.179   0.8014    
# 3a - 6r == 0        0.076301   0.116031   0.658   1.0000    
# 3b - 6r == 0        0.314036   0.100352   3.129   0.1661    
# 1 - 6r == 0         0.222098   0.090035   2.467   0.5958    
# 43 - 6r == 0        0.177619   0.107017   1.660   0.9817    
# OP4 - 6r == 0       0.085794   0.095465   0.899   1.0000    
# PFop - 6r == 0      0.212076   0.096847   2.190   0.7924    
# PF - 6r == 0        0.225411   0.103934   2.169   0.8084    
# PSL - 6r == 0       0.266931   0.123039   2.169   0.8078    
# STV - 6r == 0       0.264276   0.116031   2.278   0.7360    
# A5 - 6r == 0        0.296731   0.098456   3.014   0.2210    
# A4 - 6r == 0        0.232669   0.091853   2.533   0.5445    
# 55b - 6v == 0      -0.022424   0.076322  -0.294   1.0000    
# 4 - 6v == 0        -0.195416   0.091068  -2.146   0.8191    
# 3a - 6v == 0       -0.356245   0.099544  -3.579   0.0457 *  
#   3b - 6v == 0       -0.118510   0.080723  -1.468   0.9955    
# 1 - 6v == 0        -0.210448   0.067468  -3.119   0.1696    
# 43 - 6v == 0       -0.254927   0.088873  -2.868   0.3007    
# OP4 - 6v == 0      -0.346752   0.074561  -4.651    <0.01 ***
#   PFop - 6v == 0     -0.220470   0.076322  -2.889   0.2900    
# PF - 6v == 0       -0.207135   0.085135  -2.433   0.6217    
# PSL - 6v == 0      -0.165615   0.107631  -1.539   0.9921    
# STV - 6v == 0      -0.168269   0.099544  -1.690   0.9773    
# A5 - 6v == 0       -0.135815   0.078353  -1.733   0.9712    
# A4 - 6v == 0       -0.199877   0.069876  -2.860   0.3077    
# 4 - 55b == 0       -0.172992   0.096847  -1.786   0.9609    
# 3a - 55b == 0      -0.333822   0.104857  -3.184   0.1438    
# 3b - 55b == 0      -0.096086   0.087191  -1.102   0.9999    
# 1 - 55b == 0       -0.188024   0.075087  -2.504   0.5651    
# 43 - 55b == 0      -0.232504   0.094786  -2.453   0.6050    
# OP4 - 55b == 0     -0.324328   0.081519  -3.979   0.0119 *  
#   PFop - 55b == 0    -0.198046   0.083133  -2.382   0.6601    
# PF - 55b == 0      -0.184711   0.091291  -2.023   0.8828    
# PSL - 55b == 0     -0.143192   0.112563  -1.272   0.9993    
# STV - 55b == 0     -0.145846   0.104857  -1.391   0.9977    
# A5 - 55b == 0      -0.113391   0.085002  -1.334   0.9987    
# A4 - 55b == 0      -0.177454   0.077257  -2.297   0.7218    
# 3a - 4 == 0        -0.160830   0.116031  -1.386   0.9978    
# 3b - 4 == 0         0.076906   0.100352   0.766   1.0000    
# 1 - 4 == 0         -0.015032   0.090035  -0.167   1.0000    
# 43 - 4 == 0        -0.059512   0.107017  -0.556   1.0000    
# OP4 - 4 == 0       -0.151337   0.095465  -1.585   0.9887    
# PFop - 4 == 0      -0.025055   0.096847  -0.259   1.0000    
# PF - 4 == 0        -0.011719   0.103934  -0.113   1.0000    
# PSL - 4 == 0        0.029800   0.123039   0.242   1.0000    
# STV - 4 == 0        0.027146   0.116031   0.234   1.0000    
# A5 - 4 == 0         0.059601   0.098456   0.605   1.0000    
# A4 - 4 == 0        -0.004462   0.091853  -0.049   1.0000    
# 3b - 3a == 0        0.237736   0.108102   2.199   0.7879    
# 1 - 3a == 0         0.145797   0.098600   1.479   0.9951    
# 43 - 3a == 0        0.101318   0.114317   0.886   1.0000    
# OP4 - 3a == 0       0.009493   0.103582   0.092   1.0000    
# PFop - 3a == 0      0.135775   0.104857   1.295   0.9991    
# PF - 3a == 0        0.149111   0.111435   1.338   0.9986    
# PSL - 3a == 0       0.190630   0.129438   1.473   0.9953    
# STV - 3a == 0       0.187976   0.122796   1.531   0.9926    
# A5 - 3a == 0        0.220431   0.106344   2.073   0.8596    
# A4 - 3a == 0        0.156368   0.100262   1.560   0.9908    
# 1 - 3b == 0        -0.091938   0.079556  -1.156   0.9998    
# 43 - 3b == 0       -0.136418   0.098364  -1.387   0.9977    
# OP4 - 3b == 0      -0.228243   0.085653  -2.665   0.4437    
# PFop - 3b == 0     -0.101961   0.087191  -1.169   0.9998    
# PF - 3b == 0       -0.088625   0.095001  -0.933   1.0000    
# PSL - 3b == 0      -0.047106   0.115592  -0.408   1.0000    
# STV - 3b == 0      -0.049760   0.108102  -0.460   1.0000    
# A5 - 3b == 0       -0.017305   0.088974  -0.194   1.0000    
# A4 - 3b == 0       -0.081368   0.081608  -0.997   1.0000    
# 43 - 1 == 0        -0.044479   0.087814  -0.507   1.0000    
# OP4 - 1 == 0       -0.136304   0.073295  -1.860   0.9424    
# PFop - 1 == 0      -0.010022   0.075087  -0.133   1.0000    
# PF - 1 == 0         0.003313   0.084029   0.039   1.0000    
# PSL - 1 == 0        0.044833   0.106758   0.420   1.0000    
# STV - 1 == 0        0.042179   0.098600   0.428   1.0000    
# A5 - 1 == 0         0.074633   0.077150   0.967   1.0000    
# A4 - 1 == 0         0.010571   0.068524   0.154   1.0000    
# OP4 - 43 == 0      -0.091825   0.093374  -0.983   1.0000    
# PFop - 43 == 0      0.034457   0.094786   0.364   1.0000    
# PF - 43 == 0        0.047793   0.102016   0.468   1.0000    
# PSL - 43 == 0       0.089312   0.121424   0.736   1.0000    
# STV - 43 == 0       0.086658   0.114317   0.758   1.0000    
# A5 - 43 == 0        0.119113   0.096429   1.235   0.9995    
# A4 - 43 == 0        0.055050   0.089677   0.614   1.0000    
# PFop - OP4 == 0     0.126282   0.081519   1.549   0.9913    
# PF - OP4 == 0       0.139618   0.089823   1.554   0.9910    
# PSL - OP4 == 0      0.181137   0.111376   1.626   0.9853    
# STV - OP4 == 0      0.178483   0.103582   1.723   0.9723    
# A5 - OP4 == 0       0.210938   0.083423   2.529   0.5454    
# A4 - OP4 == 0       0.146875   0.075517   1.945   0.9146    
# PF - PFop == 0      0.013336   0.091291   0.146   1.0000    
# PSL - PFop == 0     0.054855   0.112563   0.487   1.0000    
# STV - PFop == 0     0.052201   0.104857   0.498   1.0000    
# A5 - PFop == 0      0.084655   0.085002   0.996   1.0000    
# A4 - PFop == 0      0.020593   0.077257   0.267   1.0000    
# PSL - PF == 0       0.041519   0.118715   0.350   1.0000    
# STV - PF == 0       0.038865   0.111435   0.349   1.0000    
# A5 - PF == 0        0.071320   0.092995   0.767   1.0000    
# A4 - PF == 0        0.007257   0.085974   0.084   1.0000    
# STV - PSL == 0     -0.002654   0.129438  -0.021   1.0000    
# A5 - PSL == 0       0.029801   0.113950   0.262   1.0000    
# A4 - PSL == 0      -0.034262   0.108296  -0.316   1.0000    
# A5 - STV == 0       0.032455   0.106344   0.305   1.0000    
# A4 - STV == 0      -0.031608   0.100262  -0.315   1.0000    
# A4 - A5 == 0       -0.064063   0.079264  -0.808   1.0000    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)


