library(tidyverse)
library(ggExtra)
library(coin)

theme_set(theme_bw())

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')

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

# DBS_aper <- read_csv('data/DBS_aperiodic_components_all_avg.csv') %>%
#   select(-X1) %>%
#   mutate(session_id=ceiling(session_id)) %>%
#   mutate(electrode=tolower(electrode)) %>%
#   left_join(subjects %>% select(subject,dx,dbs_target)) %>%
#   filter(!subject%in%c('DBS3014')) #strange noise with peaks for DBS recordings

DBS_aper_coord_dx <- read_tsv('data/A01_DBS_aper_coord_dx.tsv') %>%
  mutate(electrode=tolower(electrode)) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel) %>%
  filter(!str_detect(electrode,fixed('dbs_r'))) 


color_pd_ecog = '#6F67A6' #PD ECoG
color_pd_stn = '#F7924A' #PD STN
#=====

f_min = 1

aper_curves <- DBS_aper_coord_dx %>%
  filter(dbs_target=='STN') %>%
  filter(dx%in%c('PD')) %>%
  mutate(id=paste(subject,session_id,electrode,sep='_')) %>%
  group_by(id,subject,electrode,session_id,electrode_type,dx,dbs_target) %>%
  group_modify(~{ 
    tibble(f=10^seq(0,2.4,0.1)) %>%
      mutate(P= .x$offset * ((10^(.x$log_knee*.x$exponent) + f_min^.x$exponent)/(10^(.x$log_knee*.x$exponent) + f^.x$exponent)))
  })

aper_curves %>%
  ggplot() +
  aes(x=f,y=P,color=electrode_type,alpha=electrode_type,group=id)+
  geom_line(size=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values=c(color_pd_stn,color_pd_ecog),labels=c('STN','Cortex')) + 
  scale_alpha_manual(values=c(0.12,0.03),labels=c('STN','Cortex'))+
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position=c(0.2,0.2), legend.background = element_blank()) +
  coord_cartesian(ylim=c(0.0001,10),xlim=c(1,200),expand=FALSE) +
  labs(color='Region',alpha='Region')
ggsave('fig/B07_01_aper_curve_STN-vs-cortex.pdf',width=4,height=4)

#=============
DBS_aper_agg <- DBS_aper_coord_dx %>%
  filter(electrode_type=='dbs') %>%
  filter(dbs_target=='STN') %>%
  group_by(subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) %>%
  select(subject,offset,log_knee,exponent,updrs_preop_on,updrs_preop_off) %>%
  pivot_longer(cols=c(updrs_preop_on,updrs_preop_off),names_to='UPDRS',values_to='updrs_score') %>%
  mutate(UPDRS=str_remove(UPDRS,fixed('updrs_preop_')))

palette <- c('#80666C','#EB8706')

DBS_aper_agg %>%
  ggplot() +
  aes(x=updrs_score,y=exponent,color=UPDRS)+
  geom_point() +
  scale_color_manual(values=palette)
ggsave('fig/B07_02_exponent-vs-UPDRS.pdf',width=4,height=3)

DBS_aper_agg %>%
  ggplot() +
  aes(x=updrs_score,y=log_knee,color=UPDRS)+
  geom_point() 
  
DBS_aper_agg %>%
  ggplot() +
  aes(x=updrs_score,y=offset,color=UPDRS)+
  geom_point() 

#-------------

DBS_aper_agg %>% 
  filter(UPDRS=='on') %>%
  with(cor.test(exponent,updrs_score,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  exponent and updrs_score
# S = 159.98, p-value = 0.9336
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#        rho 
# 0.03039528 

DBS_aper_agg %>% 
  filter(UPDRS=='off') %>%
  with(cor.test(exponent,updrs_score,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  exponent and updrs_score
# S = 238, p-value = 0.6037
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1678322 

DBS_aper_agg %>%
  ggplot() +
  aes(x=updrs_score,y=offset,color=UPDRS)+
  geom_point() +
  scale_color_manual(values=palette)
ggsave('fig/B07_02b_offset-vs-UPDRS.pdf',width=4,height=3)

#-------------

DBS_aper_agg %>% 
  filter(UPDRS=='on') %>%
  with(cor.test(offset,updrs_score,method = 'pearson'))
# Pearson's product-moment correlation
# 
# data:  offset and updrs_score
# t = -2.1023, df = 8, p-value = 0.06869
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.89136869  0.05299566
# sample estimates:
#        cor 
# -0.5965355 


DBS_aper_agg %>% 
  filter(UPDRS=='off') %>%
  with(cor.test(offset,updrs_score,method = 'pearson'))
# Pearson's product-moment correlation
# 
# data:  offset and updrs_score
# t = -0.26702, df = 10, p-value = 0.7949
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.6277296  0.5146114
# sample estimates:
#         cor 
# -0.08413977 
#=============



DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog' & dbs_target=='GPi') %>%
  with(unique(subject))

DBS_aper %>%
  filter(electrode_type=='ecog' & dbs_target=='GPi') %>%
  with(unique(subject))



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

with(DBS_aper_agg2,unique(electrode_type_dbs_target))

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG
color_ep_seeg = '#8A4F80' #EP sEEG
color_pd_stn = '#F7924A' #PD STN
color_pd_gpi = '#F9BD00' #PD GPi
color_et_vim = '#36A5D1' #ET VIM
color_ep_cm = '#9EB859' #EP CM
palette = c(color_et_vim,color_et_ecog,color_pd_stn,color_pd_ecog,color_pd_gpi,color_pd_ecog)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,ymin=exponent-sd_exponent,ymax=exponent+sd_exponent,x=electrode_type,color=electrode_type_dbs_target)+
  #geom_errorbar(position=position_dodge2(width = 0.5)) +
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  facet_wrap(~dbs_target) + 
  scale_color_manual(values=palette) +
  scale_x_discrete(labels=c("DBS","ECoG")) +
  labs(y='Exponent')
ggsave('fig/B07_03_exponent-vs-electrode-type.pdf',width=4,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=log_knee,ymin=log_knee-sd_log_knee,ymax=log_knee+sd_log_knee,x=electrode_type,color=electrode_type_dbs_target)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  geom_hline(yintercept=0,color='gray',linetype=2) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  facet_wrap(~dbs_target) + 
  scale_color_manual(values=palette) +
  scale_x_discrete(labels=c("DBS","ECoG")) +
  coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Knee Frequency (Hz)')
ggsave('fig/B07_04_log-knee-vs-electrode-type.pdf',width=4,height=3)


DBS_aper_agg2 %>%
  ggplot() +
  aes(y=offset,ymin=offset-sd_offset,ymax=offset+sd_offset,x=electrode_type,color=electrode_type_dbs_target)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  #scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  facet_wrap(~dbs_target) + 
  scale_color_manual(values=palette) +
  scale_x_discrete(labels=c("DBS","ECoG")) +
  #coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Offset')
ggsave('fig/B07_05_offset-vs-electrode-type.pdf',width=4,height=3)





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


within_test_db <- DBS_aper_coord_dx %>%
  filter(tolower(type)%in%c('ecog','dbs'))%>%
  filter(toupper(dbs_target)%in%c('STN','GPI','VIM')) %>%
  filter(toupper(dx)%in%c('PD','ET')) %>% 
  group_by(subject) %>%
  group_modify(oneway_test_pvalue) %>%
  pivot_longer(cols=c(exponent,offset,log_knee),values_to='pvalue',names_to = 'param') %>%
  ungroup() %>%
  mutate(pvalue_fdr = p.adjust(pvalue,method="fdr")) %>% 
  mutate(is_signif = pvalue_fdr < 0.05) %>%
  write_tsv('fig/B07_03_pvalues-fdr.tsv') %>%
  pivot_wider(id_cols=c(subject),names_from=param,values_from = is_signif, values_fill=FALSE, names_prefix = 'is_signif_')

DBS_aper_agg2 %>%
  left_join(within_test_db) %>%
  ggplot() +
  aes(y=exponent,ymin=exponent-sd_exponent,ymax=exponent+sd_exponent,
      x=electrode_type,color=electrode_type_dbs_target,linetype=is_signif_exponent)+
  #geom_errorbar(position=position_dodge2(width = 0.5)) +
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(aes(color=electrode_type_dbs_target,linetype=NULL),position=position_dodge2(width = 0.75),size=0.2) +
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  facet_wrap(~dbs_target) + 
  scale_linetype_discrete() +
  scale_color_manual(values=palette) +
  scale_x_discrete(labels=c("DBS","ECoG")) +
  labs(y='Exponent')
ggsave('fig/B07_03b_exponent-vs-electrode-type.pdf',width=4,height=3)

DBS_aper_agg2 %>%
  left_join(within_test_db) %>%
  ggplot() +
  aes(y=log_knee,ymin=log_knee-sd_log_knee,ymax=log_knee+sd_log_knee,x=electrode_type,color=electrode_type_dbs_target)+
  geom_line(aes(group = subject, linetype=is_signif_log_knee),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  geom_hline(yintercept=0,color='gray',linetype=2) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  facet_wrap(~dbs_target) + 
  scale_color_manual(values=palette) +
  scale_x_discrete(labels=c("DBS","ECoG")) +
  coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Knee Frequency (Hz)')
ggsave('fig/B07_04b_log-knee-vs-electrode-type.pdf',width=4,height=3)


DBS_aper_agg2 %>%
  left_join(within_test_db) %>%
  ggplot() +
  aes(y=offset,ymin=offset-sd_offset,ymax=offset+sd_offset,x=electrode_type,color=electrode_type_dbs_target)+
  geom_line(aes(group = subject, linetype=is_signif_offset),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  #scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  facet_wrap(~dbs_target) + 
  scale_color_manual(values=palette) +
  scale_x_discrete(labels=c("DBS","ECoG")) +
  #coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Offset')
ggsave('fig/B07_05b_offset-vs-electrode-type.pdf',width=4,height=3)



#---------
#STN

DBS_aper_agg2_for_test <- DBS_aper_agg2 %>%
  filter(dbs_target=='STN') %>%
  group_by(subject) %>%
  mutate(has_dbs='dbs'%in%electrode_type) %>%
  filter(has_dbs) %>%
  select(subject, electrode_type,offset,log_knee,exponent) %>%
  pivot_wider(names_from=electrode_type, values_from=c(offset,log_knee,exponent))

DBS_aper_agg2_for_test %>%
  with(wilcox.test(offset_dbs,offset_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  offset_dbs and offset_ecog
# V = 9, p-value = 0.008057
# alternative hypothesis: true location shift is not equal to 0

DBS_aper_agg2_for_test %>%
  with(wilcox.test(exponent_dbs,exponent_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  exponent_dbs and exponent_ecog
# V = 0, p-value = 0.0002441
# alternative hypothesis: true location shift is not equal to 0

DBS_aper_agg2_for_test %>%
  with(wilcox.test(log_knee_dbs,log_knee_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  log_knee_dbs and log_knee_ecog
# V = 0, p-value = 0.0002441
# alternative hypothesis: true location shift is not equal to 0


#---
# VIM

DBS_aper_agg2_for_test <- DBS_aper_agg2 %>%
  filter(dbs_target=='VIM') %>%
  group_by(subject) %>%
  mutate(has_dbs='dbs'%in%electrode_type) %>%
  filter(has_dbs) %>%
  select(subject, electrode_type,offset,log_knee,exponent) %>%
  pivot_wider(names_from=electrode_type, values_from=c(offset,log_knee,exponent))


DBS_aper_agg2_for_test %>%
  with(wilcox.test(offset_dbs,offset_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  offset_dbs and offset_ecog
# V = 16, p-value = 0.01025
# alternative hypothesis: true location shift is not equal to 0

DBS_aper_agg2_for_test %>%
  with(wilcox.test(exponent_dbs,exponent_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  exponent_dbs and exponent_ecog
# V = 0, p-value = 6.104e-05
# alternative hypothesis: true location shift is not equal to 0

DBS_aper_agg2_for_test %>%
  with(wilcox.test(log_knee_dbs,log_knee_ecog, paired=TRUE))

# Wilcoxon signed rank exact test
# 
# data:  log_knee_dbs and log_knee_ecog
# V = 0, p-value = 6.104e-05
# alternative hypothesis: true location shift is not equal to 0
# 


#---
#GPi

DBS_aper_agg2_for_test <- DBS_aper_agg2 %>%
  filter(dbs_target=='GPi') %>%
  group_by(subject) %>%
  mutate(has_dbs='dbs'%in%electrode_type) %>%
  filter(has_dbs) %>%
  select(subject, electrode_type,offset,log_knee,exponent) %>%
  pivot_wider(names_from=electrode_type, values_from=c(offset,log_knee,exponent))



DBS_aper_agg2_for_test %>%
  with(wilcox.test(offset_dbs,offset_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  offset_dbs and offset_ecog
# V = 0, p-value = 0.25
# alternative hypothesis: true location shift is not equal to 0

DBS_aper_agg2_for_test %>%
  with(wilcox.test(exponent_dbs,exponent_ecog, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  exponent_dbs and exponent_ecog
# V = 0, p-value = 0.25
# alternative hypothesis: true location shift is not equal to 0


DBS_aper_agg2_for_test %>%
  with(wilcox.test(log_knee_dbs,log_knee_ecog, paired=TRUE))

# Wilcoxon signed rank exact test
# 
# data:  log_knee_dbs and log_knee_ecog
# V = 0, p-value = 0.25
# alternative hypothesis: true location shift is not equal to 0

####



DBS_aper_coord_dx %>%
  filter(tolower(type)%in%c('ecog','dbs'))%>%
  filter(toupper(dbs_target)%in%c('STN','GPI','VIM')) %>%
  filter(toupper(dx)%in%c('PD','ET')) %>%
  group_by(subject) %>%
  mutate(has_dbs='dbs'%in%electrode_type) %>%
  filter(has_dbs) %>%
  group_by(dbs_target,electrode_type) %>%
  summarise(N_knee=sum(log_knee>0),N_no_knee=sum(log_knee<0),N=n())

# dbs_target electrode_type N_knee N_no_knee     N
# <chr>      <chr>           <int>     <int> <int>
# 1 GPi        dbs                 6        16    22
# 2 GPi        ecog              480         8   488
# 3 STN        dbs                11        62    73
# 4 STN        ecog             3835        42  3877
# 5 VIM        dbs                22       151   173
# 6 VIM        ecog             2263        20  2283

  group_by(dbs_target) %>%
  
chisq.test(matrix(c(6,420,16,8),nrow=2))
  # Pearson's Chi-squared test with Yates' continuity correction
  # 
  # data:  matrix(c(6, 420, 16, 8), nrow = 2)
  # X-squared = 194.29, df = 1, p-value < 2.2e-16 
  fisher.test(matrix(c(6,420,16,8),nrow=2))
#   Fisher's Exact Test for Count Data
# 
# data:  matrix(c(6, 420, 16, 8), nrow = 2)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.001836983 0.026222544
# sample estimates:
#  odds ratio 
# 0.007591096 
  
  # 
  select(subject, electrode_type,offset,log_knee,exponent) %>%
  pivot_wider(names_from=electrode_type, values_from=c(offset,log_knee,exponent))



DBS_aper_agg2_for_test %>%
  with(table(factor(log_knee_dbs > 0),factor(log_knee_ecog > 0)))



with(fisher.test(log_knee_dbs > 0,log_knee_ecog > 0))

table



summary(lm(exponent ~ electrode_type + dbs_target, data=DBS_aper_agg2))


oneway_test(exponent ~ factor(electrode_type), data=DBS_aper_agg2, distribution=approximate(nresample=999))


summary(lm(log_knee ~ electrode_type, data=DBS_aper_agg2))


oneway_test(log_knee ~ factor(electrode_type), data=DBS_aper_agg2, distribution=approximate(nresample=999))


summary(lm(offset ~ electrode_type, data=DBS_aper_agg2))

#=======================

DBS_aper_summary_stats <- DBS_aper_agg2 %>%
  group_by(dbs_target, electrode_type) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
          sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
          offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  write_tsv('data/DBS_aper_summary_stats.tsv')


#=======================

DBS_aper %>%
  filter(electrode_type%in%c('ecog','dbs'))%>%
  filter(dbs_target=='STN') %>%
  filter(dx%in%c('PD')) %>%
  ggplot() +
  aes(y=exponent,x=log_knee,color=electrode_type) +
  geom_density2d(alpha=0.5,breaks=c(0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.3)+
  geom_vline(xintercept=0,color='gray',linetype=2,size=0.8)+
  theme(legend.position=c(0.2,0.8),legend.background = element_blank())+
  labs(color="Region",x="Knee Frequency (Hz)") +
  scale_x_continuous(limits=c(-1,1.7),breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.3,1,3,10,30))+
  scale_y_continuous(breaks=1:5,limits=c(0.7,5)) +
  scale_color_manual(values=c(color_pd_stn,color_pd_ecog),labels=c('STN','Cortex'))+
  labs(y='Exponenet')
ggsave('fig/B07_06_exponent-vs-log-knee_STN-vs-cortex.pdf',width=4,height=4)

#===========
