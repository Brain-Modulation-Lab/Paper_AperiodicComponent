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

color_ep_cm = '#36A5D1' #EP CM
color_ep_seeg = '#C4604F' #EP sEEG

#=====

f_min = 1

aper_curves <- sEEG_aper %>%
  mutate(id=paste(subject,session_id,electrode,sep='_')) %>%
  group_by(id,subject,electrode,session_id,localization) %>%
  group_modify(~{ 
    tibble(f=10^seq(0,2.4,0.1)) %>%
      mutate(P= .x$offset * ((10^(.x$log_knee*.x$exponent) + f_min^.x$exponent)/(10^(.x$log_knee*.x$exponent) + f^.x$exponent)))
  })

aper_curves %>%
  ggplot() +
  aes(x=f,y=P,color=localization,alpha=localization,group=id)+
  geom_line(size=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values=c(color_ep_seeg,color_ep_cm),labels=c('Cortex','CM')) + 
  scale_alpha_manual(values=c(0.15,0.25),labels=c('Cortex','CM'))+
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position=c(0.2,0.2), legend.background = element_blank()) +
  coord_cartesian(xlim=c(1,200),expand=FALSE) +
  labs(color='Region',alpha='Region') 
ggsave('fig/B05_01_aper_curve_CM-vs-cortex.pdf',width=4,height=4)


#=============

se <- function(x) sd(x)/sqrt(length(x))

sEEG_aper_agg2 <- sEEG_aper %>%
  group_by(subject,electrode_type) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
            sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
            offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) 

sEEG_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,ymin=exponent-sd_exponent,ymax=exponent+sd_exponent,x=electrode_type,color=electrode_type)+
  #geom_errorbar(position=position_dodge2(width = 0.5)) +
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  #scale_y_continuous(limits=c(2.5,4))+
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  scale_x_discrete(labels=c("Thal","Ctx")) +
  labs(y='Exponent')
ggsave('fig/B05_02_exponent-vs-electrode-type_EP.pdf',width=2,height=3)

sEEG_aper_agg2 %>%
  ggplot() +
  aes(y=log_knee,ymin=log_knee-sd_log_knee,ymax=log_knee+sd_log_knee,x=electrode_type,color=electrode_type)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  geom_hline(yintercept=0,color='gray',linetype=2) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  scale_x_discrete(labels=c("Thal","Ctx")) +
  coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Knee Frequency (Hz)')
ggsave('fig/B05_03_log-knee-vs-electrode-type_EP.pdf',width=2.2,height=3)


sEEG_aper_agg2 %>%
  ggplot() +
  aes(y=offset,ymin=offset-sd_offset,ymax=offset+sd_offset,x=electrode_type,color=electrode_type)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  #scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  scale_x_discrete(labels=c("Thal","Ctx")) +
  #coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Offset')
ggsave('fig/B05_04_offset-vs-electrode-type_EP.pdf',width=2,height=3)

#======================
### individual subject level stats

oneway_test_pvalue<- function(d,y){
  if (length(unique(d$electrode_type))<2){
    return(tibble())  
  } else {
    e<-oneway_test(exponent ~ factor(electrode_type), distribution=approximate(nresample=9999),data=d)
    o<-oneway_test(offset ~ factor(electrode_type), distribution=approximate(nresample=9999),data=d)
    k<-oneway_test(log_knee ~ factor(electrode_type), distribution=approximate(nresample=9999),data=d)
    e<-max(as.numeric(pvalue(e)),1e-4)
    o<-max(as.numeric(pvalue(o)),1e-4)
    k<-max(as.numeric(pvalue(k)),1e-4)
    return(tibble(exponent=e,offset=o,log_knee=k))
  }
}


within_test_db <- sEEG_aper %>%
  group_by(subject) %>%
  group_modify(oneway_test_pvalue) %>%
  pivot_longer(cols=c(exponent,offset,log_knee),values_to='pvalue',names_to = 'param') %>%
  ungroup() %>%
  mutate(pvalue_fdr = p.adjust(pvalue,method="fdr")) %>% 
  mutate(is_signif = pvalue_fdr < 0.05) %>%
  write_tsv('fig/B05_02_pvalues-fdr.tsv') %>%
  pivot_wider(id_cols=c(subject),names_from=param,values_from = is_signif, values_fill=FALSE, names_prefix = 'is_signif_')


sEEG_aper_agg2 %>%
  left_join(within_test_db) %>%
  ggplot() +
  aes(y=exponent,ymin=exponent-sd_exponent,ymax=exponent+sd_exponent,x=electrode_type,color=electrode_type)+
  #geom_errorbar(position=position_dodge2(width = 0.5)) +
  geom_line(aes(group = subject,linetype=is_signif_exponent),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  #scale_y_continuous(limits=c(2.5,4))+
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  scale_x_discrete(labels=c("Thal","Ctx")) +
  labs(y='Exponent')
ggsave('fig/B05_02b_exponent-vs-electrode-type_EP.pdf',width=2,height=3)

sEEG_aper_agg2 %>%
  left_join(within_test_db) %>%
  ggplot() +
  aes(y=log_knee,ymin=log_knee-sd_log_knee,ymax=log_knee+sd_log_knee,x=electrode_type,color=electrode_type)+
  geom_line(aes(group = subject,linetype=is_signif_log_knee),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  geom_hline(yintercept=0,color='gray',linetype=2) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  scale_x_discrete(labels=c("Thal","Ctx")) +
  coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Knee Frequency (Hz)')
ggsave('fig/B05_03b_log-knee-vs-electrode-type_EP.pdf',width=2.2,height=3)


sEEG_aper_agg2 %>%
  left_join(within_test_db) %>%
  ggplot() +
  aes(y=offset,ymin=offset-sd_offset,ymax=offset+sd_offset,x=electrode_type,color=electrode_type)+
  geom_line(aes(group = subject, linetype=is_signif_offset),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  #scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  scale_x_discrete(labels=c("Thal","Ctx")) +
  #coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Offset')
ggsave('fig/B05_04b_offset-vs-electrode-type_EP.pdf',width=2,height=3)







#=======================

sEEG_aper_summary_stats <- sEEG_aper_agg2 %>%
  group_by(dbs_target, electrode_type) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
            sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
            offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  write_tsv('data/sEEG_aper_summary_stats.tsv')


sEEG_aper_agg2_for_test <- sEEG_aper_agg2 %>%
  group_by(subject) %>%
  select(subject, electrode_type,offset,log_knee,exponent) %>%
  pivot_wider(names_from=electrode_type, values_from=c(offset,log_knee,exponent))


sEEG_aper_agg2_for_test %>%
  with(wilcox.test(offset_subcortical,offset_cortical, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  offset_subcortical and offset_cortical
# V = 0, p-value = 0.007813
# alternative hypothesis: true location shift is not equal to 0

sEEG_aper_agg2_for_test %>%
  with(wilcox.test(exponent_subcortical,exponent_cortical, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  exponent_subcortical and exponent_cortical
# V = 0, p-value = 0.007813
# alternative hypothesis: true location shift is not equal to 0

sEEG_aper_agg2_for_test %>%
  with(wilcox.test(log_knee_subcortical,log_knee_cortical, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  log_knee_subcortical and log_knee_cortical
# V = 0, p-value = 0.007813
# alternative hypothesis: true location shift is not equal to 0


#=======================

sEEG_aper %>%
  ggplot() +
  aes(y=exponent,x=log_knee,color=electrode_type) +
  geom_density2d(alpha=0.5,breaks=c(0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.3)+
  geom_vline(xintercept=0,color='gray',linetype=2,size=0.8)+
  theme(legend.position=c(0.2,0.8),legend.background = element_blank())+
  labs(color="Region",x="Knee Frequency (Hz)") +
  scale_x_continuous(limits=c(-1,1.7),breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.3,1,3,10,30))+
  scale_y_continuous(breaks=1:4,limits=c(0.7,4.3)) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg),labels=c('CM','Cortex'))+
  labs(y='Exponenet')
ggsave('fig/B05_05_exponent-vs-log-knee_CM-vs-cortex.pdf',width=4,height=4)

#===========



sEEG_aper %>%
  group_by(subject) %>%
  mutate(has_thal='subcortical'%in%localization) %>%
  filter(has_thal) %>%
  group_by(localization) %>%
  summarise(N_knee=sum(log_knee>0),N_no_knee=sum(log_knee<0),N=n()) %>%
  mutate(f_no_k=N_no_knee/N,err_percent=sqrt(f_no_k*(1-f_no_k)/N))


#   localization N_knee N_no_knee     N f_no_k err_percent
#   <chr>         <int>     <int> <int>  <dbl>       <dbl>
# 1 cortical        131         5   136 0.0368      0.0161
# 2 subcortical      26        96   122 0.787       0.0371
