library(tidyverse)
library(ggExtra)

theme_set(theme_bw())

#setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')
setwd('~/Dropbox/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')

DBS_aper_new <- read_csv('data/DBS_aperiodic_components_all_avg.csv') %>%
  select(-`...1`) %>%
  mutate(session_id=ceiling(session_id)) %>%
  mutate(electrode=tolower(electrode)) %>%
  left_join(subjects %>% select(subject,dx,dbs_target))

DBS_aper_o <- read_csv('data/DBS_old_fooof_params.csv') %>%
  #select(-`...1`,-fooof_export_path) %>%
  mutate(session_id=ceiling(session_id)) %>%
  mutate(electrode=tolower(electrode)) %>%
  mutate(log_knee_o = log10(knee)) %>%
  rename(offset_o = offset, exponent_o = exponent,  knee_o = knee,
         tau_o=tau_ms,rsquared_o=rsquared) 

DBS_aper <- DBS_aper_new %>%
  left_join(DBS_aper_o)

DBS_aper_o2 <- read_csv('data/DBS_old_fooof_params.csv') %>%
  select(-`...1`,-fooof_export_path) %>%
  mutate(session_id=ceiling(session_id)) %>%
  mutate(electrode=tolower(electrode)) %>%
  mutate(param='original') %>%
  mutate(log_knee = ifelse(knee>0,log10(abs(knee)),-2)) %>%
  left_join(subjects %>% select(subject,dx,dbs_target))

DBS_aper2 <- DBS_aper_new %>%
  mutate(param='proposed') %>%
  bind_rows(DBS_aper_o2)

#=============

fig_B01_01 <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(x=rsquared_o,y=rsquared)+
  geom_point(alpha=0.5) +
  coord_equal() +
  labs(x='Original Parameterizarion R^2',y='Proposed Parameterization R^2')
fig_B01_01 <- ggMarginal(fig_B01_01, type="histogram")
fig_B01_01
ggsave('fig/B01_01_orig_vs_new_R2_corr_marginals.pdf',fig_B01_01,width=4,height=4)

#=============

DBS_aper2 %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(y=exponent,x=offset,color=param) +
  geom_density2d(alpha=0.5,breaks=c(0.05,0.1,0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.3)+
  scale_x_continuous(breaks=0:9)+
  scale_y_continuous(limits=c(1,5),breaks=1:5) +
  theme(legend.position=c(0.8,0.2),legend.background = element_blank())+
  labs(color="Parametrization")
ggsave('fig/B01_02_orig_vs_new__exponent_vs_offset.pdf',width=4,height=4)


DBS_aper2 %>%
  filter(electrode_type=='dbs') %>%
  filter(dx=='PD') %>%
  filter(dbs_target=='STN') %>%
  ggplot() +
  aes(y=exponent,x=offset,color=param) +
  geom_density2d(alpha=0.5,breaks=c(0.05,0.1,0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.6)+
  scale_x_continuous(limits=c(-1.5,3))+
  scale_y_continuous(limits=c(0,2.5),breaks=seq(0,2.5,0.5)) +
  theme(legend.position=c(0.8,0.2),legend.background = element_blank())+
  labs(color="Parametrization")
ggsave('fig/B01_02_orig_vs_new__exponent_vs_offset_DBS.pdf',width=4,height=4)


DBS_aper2 %>%
  filter(electrode_type=='dbs') %>%
  filter(param=="proposed") %>%
  filter(dx=='PD') %>%
  filter(dbs_target=='STN') %>%
  ggplot() +
  aes(y=exponent,x=log_knee,color=param) +
  #geom_density2d(alpha=0.5,breaks=c(0.05,0.1,0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.6)+
  #scale_x_continuous(limits=c(-1.5,3))+
  geom_smooth(method='lm') +
  #scale_y_continuous(limits=c(0,2.5),breaks=seq(0,2.5,0.5)) +
  theme(legend.position=c(0.8,0.2),legend.background = element_blank())+
  labs(color="Parametrization")
ggsave('fig/B01_02_orig_vs_new__exponent_vs_knee_DBS.pdf',width=4,height=4)


DBS_aper2 %>%
  filter(electrode_type=='dbs') %>%
  filter(param=="proposed") %>%
  filter(dx=='PD') %>%
  filter(dbs_target=='STN') %>%
  filter(subject!='DBS3014') %>%
  with(cor.test(exponent,log_knee,method="spearman"))

# Spearman's rank correlation rho
# 
# data:  exponent and log_knee
# S = 19119, p-value = 0.192
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1838487 
       




DBS_aper2 %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(y=exponent,x=log_knee,color=param) +
  geom_density2d(alpha=0.5,breaks=c(0.05,0.1,0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.3)+
  theme(legend.position=c(0.8,0.2),legend.background = element_blank())+
  labs(color="Parametrization") +
  scale_x_continuous(limits=c(0,6),breaks=0:6)+
  scale_y_continuous(limits=c(1,5),breaks=1:5) 
ggsave('fig/B01_03_orig-vs-new_exponent-vs-log-knee.pdf',width=4,height=4)


DBS_aper2 %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(y=offset,x=log_knee,color=param) +
  geom_density2d(alpha=0.5,breaks=c(0.05,0.1,0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.3)+
  theme(legend.position=c(0.8,0.2),legend.background = element_blank())+
  labs(color="Parametrization") +
  scale_x_continuous(limits=c(0,6),breaks=0:6)+
  scale_y_continuous(breaks=0:9)
ggsave('fig/B01_04_orig-vs-new_offset-vs-log-knee.pdf',width=4,height=4)

#====

fig_B01_05 <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(x=exponent_o,y=exponent)+
  geom_point(alpha=0.5) +
  coord_equal() +
  labs(x='Original Exponent',y='Novel Exponent')
fig_B01_05 <- ggMarginal(fig_B01_05, type="histogram")
fig_B01_05
ggsave('fig/B01_05_orig_vs_new_Exponent.pdf',fig_B01_05,width=4,height=4)

fig_B01_06 <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(x=knee_o,y=log_knee)+
  geom_point(alpha=0.5) +
  scale_color_manual(values=c('gray20','red')) +
  geom_hline(yintercept=0,color='gray',linetype=2)+
  labs(x='Original Knee parameter',y='Novel Knee frequency')+
  theme(legend.position = 'none')+
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  scale_x_log10()
fig_B01_06 <- ggMarginal(fig_B01_06, type="histogram")
fig_B01_06
ggsave('fig/B01_06_orig_vs_new_log-knee.pdf',fig_B01_06,width=4,height=4.2)


fig_B01_07 <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  ggplot() +
  aes(x=offset_o,y=offset)+
  geom_point(alpha=0.5) +
  labs(x='Original Offset',y='Novel Offset')+
  scale_x_continuous(limits=c(1,8))
fig_B01_07 <- ggMarginal(fig_B01_07, type="histogram")
fig_B01_07
ggsave('fig/B01_07_orig_vs_new_offset.pdf',fig_B01_07,width=4,height=4)

