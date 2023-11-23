library(coin)
library(multcomp)
library(lmerTest)
library(ggExtra)
library(tidyverse)

theme_set(theme_bw())

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')


subjects %>%
  group_by(dx) %>%
  summarise(n=n())


artifacts <- read_tsv('data/artifacts.txt') %>%
  fill(subject,session_id,type) %>%
  mutate(subject=str_replace(subject,fixed('-'),'')) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_'))

# artifacts %>% write_tsv('data/artifacts.txt')

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
  mutate(electrode=tolower(electrode)) %>%
  left_join(subjects %>% select(subject,dx,dbs_target)) %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel)

DBS_aper_coord_dx <- read_tsv('data/A01_DBS_aper_coord_dx.tsv') %>%
  mutate(subsesel = paste(subject,session_id,electrode,sep='_')) %>%
  filter(!subsesel %in% artifacts$subsesel)

MMP1_regions <- read_tsv('data/MMP1-regions.tsv')

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG

#=============
DBS_aper_agg <- DBS_aper %>%
  filter(rsquared > 0.98) %>%
  filter(!subject %in% c('DBS4080','DBS4086')) %>% #strange PSD plots
  filter(!subject %in% c('DBS4088')) %>% #epilepsy patient
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  #filter(dbs_target=='STN') %>%
  group_by(subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) %>%
  select(subject,offset,log_knee,exponent,updrs_preop_on,updrs_preop_off) %>%
  pivot_longer(cols=c(updrs_preop_on,updrs_preop_off),names_to='UPDRS',values_to='updrs_score') %>%
  mutate(UPDRS=str_remove(UPDRS,fixed('updrs_preop_')))


DBS_aper_agg %>%
  filter(UPDRS=='on') %>%
  ggplot() +
  aes(x=updrs_score,y=exponent)+
  geom_point() +
  geom_smooth(method='lm') 
ggsave('fig/B02_01a_exponent-vs-UPDRS-on.pdf',width=4,height=3)

DBS_aper_agg %>%
  filter(UPDRS=='off') %>%
  ggplot() +
  aes(x=updrs_score,y=exponent)+
  geom_point() +
  geom_smooth(method='lm') 
ggsave('fig/B02_01b_exponent-vs-UPDRS-off.pdf',width=4,height=3)

DBS_aper_agg %>%
  pivot_wider(names_from=UPDRS,values_from=updrs_score) %>% 
  mutate(on_off_percent = 100*on/off) %>%
  ggplot() +
  aes(x=on_off_percent,y=exponent)+
  geom_point() +
  geom_smooth(method='lm') 


DBS_aper_coord_dx %>%
  group_by(subject,dx) %>%
  summarise(n=n()) %>%
  group_by(dx) %>%
  summarise(ns=n())

DBS_aper_agg_sensorimotor <- DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
  filter(dbs_target=='STN') %>%
  filter(rsquared > 0.98) %>%
  filter(HCPMMP1_label_1 %in% c('1', '3b', '4','3a','6a','6d','FEF','PEF','55b','6r','6v')) %>%
  group_by(subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent),n=n()) %>%
  left_join(subjects) %>%
  select(subject,offset,log_knee,exponent,updrs_preop_on,updrs_preop_off,n) %>%
  pivot_longer(cols=c(updrs_preop_on,updrs_preop_off),names_to='UPDRS',values_to='updrs_score') %>%
  mutate(UPDRS=str_remove(UPDRS,fixed('updrs_preop_')))

write_tsv(DBS_aper_agg_sensorimotor,'data/dbs_aper_sensorimotorcortex.tsv')


DBS_aper_agg_sensorimotor %>%
  filter(UPDRS=='on') %>%
  ggplot() +
  aes(x=updrs_score,y=exponent)+
  geom_point() +
  geom_smooth(method='lm') 
ggsave('fig/B02_01c_exponent-sensorimotor-premotor-vs-UPDRS-on.pdf',width=4,height=3)

DBS_aper_agg_sensorimotor %>%
  filter(UPDRS=='off') %>%
  ggplot() +
  aes(x=updrs_score,y=exponent)+
  geom_point() +
  geom_smooth(method='lm') 
ggsave('fig/B02_01d_exponent-sensorimotor-premotor-vs-UPDRS-off.pdf',width=4,height=3)

DBS_aper_agg_sensorimotor %>%
  pivot_wider(names_from=UPDRS,values_from=updrs_score) %>% 
  mutate(on_off_percent = 100*on/off) %>%
  ggplot() +
  aes(x=on_off_percent,y=exponent)+
  geom_point() +
  geom_smooth(method='lm') 
ggsave('fig/B02_01e_exponent-sensorimotor-premotor-vs-UPDRS--on-off-percent.pdf',width=4,height=3)



#------

summary(lm(exponent~updrs_score,data=DBS_aper_agg_sensorimotor %>% filter(UPDRS=='on')))

DBS_aper_agg_sensorimotor %>% 
  filter(UPDRS=='on') %>%
  with(cor.test(exponent,updrs_score,method = 'spearman'))
# Spearman's rank correlation rho
# 
# data:  exponent and updrs_score
# S = 1213.6, p-value = 0.05833
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.4003965 

DBS_aper_agg_sensorimotor %>% 
  filter(UPDRS=='on') %>%
  with(cor.test(exponent,updrs_score,method = 'pearson'))

#Pearson's product-moment correlation
#
# data:  exponent and updrs_score
# t = 2.2436, df = 21, p-value = 0.03577
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.03361792 0.72120574
# sample estimates:
#       cor 
# 0.4397266 


summary(lm(exponent~updrs_score,data=DBS_aper_agg_sensorimotor %>% filter(UPDRS=='off')))

DBS_aper_agg_sensorimotor %>% 
  filter(UPDRS=='off') %>%
  with(cor.test(exponent,updrs_score,method = 'spearman'))
# data:  exponent and updrs_score
# S = 1838.4, p-value = 0.1553
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# 0.2929177 

DBS_aper_agg_sensorimotor %>% 
  filter(UPDRS=='off') %>%
  with(cor.test(exponent,updrs_score,method = 'pearson'))
# Pearson's product-moment correlation
# 
# data:  exponent and updrs_score
# t = 1.2657, df = 23, p-value = 0.2183
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1556532  0.5907404
# sample estimates:
#       cor 
# 0.2551716 

DBS_aper_agg_sensorimotor %>%
  pivot_wider(names_from=UPDRS,values_from=updrs_score) %>% 
  mutate(on_off_percent = 100*on/off) %>%
  with(cor.test(exponent,on_off_percent,method = 'spearman'))

# Spearman's rank correlation rho
# 
# data:  exponent and on_off_percent
# S = 1622.9, p-value = 0.3647
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.1981715 


DBS_aper_agg_sensorimotor %>%
  pivot_wider(names_from=UPDRS,values_from=updrs_score) %>% 
  mutate(on_off_percent = 100*on/off) %>%
  filter(!is.na(on_off_percent)) %>% 
  with(cor.test(exponent,on_off_percent,method = 'pearson'))

# Pearson's product-moment correlation
# 
# data:  exponent and on_off_percent
# t = 1.5781, df = 21, p-value = 0.1295
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1000264  0.6504974
# sample estimates:
#       cor 
# 0.3256007 

#=============

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG

DBS_aper_agg2 <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  group_by(subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) %>%
  mutate(dx=factor(dx,levels=c('PD','ET'),ordered=TRUE))

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,x=dx,fill=dx)+
  geom_dotplot(binaxis = "y", stackdir = "center") +
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_y_continuous(limits=c(2.5,4.1))+
  scale_fill_manual(values=c(color_pd_ecog,color_et_ecog))
ggsave('fig/B02_02_exponent-vs-dx.pdf',width=1.5,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=log_knee,x=dx,fill=dx)+
  geom_dotplot(binaxis = "y", stackdir = "center") +
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_fill_manual(values=c(color_pd_ecog,color_et_ecog))
ggsave('fig/B02_03_log-knee-vs-dx.pdf',width=1.5,height=3)

DBS_aper_agg2 %>%
  ggplot() +
  aes(y=offset,x=dx,fill=dx)+
  geom_dotplot(binaxis = "y", stackdir = "center") +
  theme(legend.position='none',panel.grid = element_blank()) +
  scale_fill_manual(values=c(color_pd_ecog,color_et_ecog))
ggsave('fig/B02_04_offset-vs-dx.pdf',width=1.5,height=3)



#---------

summary(lm(exponent ~ dx, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.24981    0.05171  62.847   <2e-16 ***
#   dx.L        -0.07850    0.07313  -1.073    0.289    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.349 on 47 degrees of freedom
# Multiple R-squared:  0.02393,	Adjusted R-squared:  0.003162 
# F-statistic: 1.152 on 1 and 47 DF,  p-value: 0.2886

oneway_test(exponent ~ factor(dx), data=DBS_aper_agg2, distribution=approximate(nresample=999))
# Approximative Linear-by-Linear Association Test
# 
# data:  exponent by factor(dx) (PD < ET)
# Z = -1.0717, p-value = 0.2833
# alternative hypothesis: two.sided

summary(lm(log_knee ~ dx, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.21017    0.02050  59.042   <2e-16 ***
#   dx.L         0.01193    0.02899   0.412    0.683    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1383 on 47 degrees of freedom
# Multiple R-squared:  0.00359,	Adjusted R-squared:  -0.01761 
# F-statistic: 0.1694 on 1 and 47 DF,  p-value: 0.6826

oneway_test(log_knee ~ factor(dx), data=DBS_aper_agg2, distribution=approximate(nresample=999))
# Approximative Linear-by-Linear Association Test
# 
# data:  log_knee by factor(dx) (PD < ET)
# Z = 0.41514, p-value = 0.6817
# alternative hypothesis: two.sided

summary(lm(offset ~ dx, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.66148    0.07936  20.936   <2e-16 ***
#   dx.L        -0.12557    0.11223  -1.119    0.269    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5356 on 47 degrees of freedom
# Multiple R-squared:  0.02594,	Adjusted R-squared:  0.00522 
# F-statistic: 1.252 on 1 and 47 DF,  p-value: 0.2689

oneway_test(offset ~ factor(dx), data=DBS_aper_agg2, distribution=approximate(nresample=999))
# Approximative Linear-by-Linear Association Test
# 
# data:  offset by factor(dx) (PD < ET)
# Z = -1.116, p-value = 0.2813
# alternative hypothesis: two.sided

#=======================

DBS_aper_agg3 <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  group_by(subject) %>%
  summarise(offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) %>%
  left_join(subjects) %>%
  mutate(dx=factor(dx,levels=c('PD','ET'),ordered=TRUE))

DBS_aper_agg3 %>%
  ggplot() +
  aes(y=exponent,x=dbs_age)+
  geom_point()+
  geom_smooth(method='lm',alpha=0.25,color='gray50')+
  theme(legend.position='none')+
  scale_y_continuous(limits=c(2.5,4))+
  scale_x_continuous(limits=c(49,80))
ggsave('fig/B02_05_exponent-vs-age.pdf',width=3,height=3)

summary(lm(exponent~dbs_age,data=DBS_aper_agg3))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 2.866221   0.562363   5.097  5.5e-05 ***
#   dbs_age     0.006616   0.008407   0.787    0.441    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3147 on 20 degrees of freedom
# Multiple R-squared:  0.03003,	Adjusted R-squared:  -0.01846 
# F-statistic: 0.6193 on 1 and 20 DF,  p-value: 0.4405

#=========

DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  mutate(dx=factor(dx,levels=c('PD','ET'),ordered=TRUE)) %>%
  ggplot() +
  aes(y=exponent,x=log_knee,color=dx) +
  geom_density2d(alpha=0.5,breaks=c(0.2,0.4,0.8))+
  geom_point(alpha=1,size=0.3)+
  geom_vline(xintercept=0,color='gray',linetype=2,size=0.8)+
  theme(legend.position=c(0.2,0.8),legend.background = element_blank())+
  labs(color="Diagnosis",x="knee frequency (Hz)") +
  scale_x_continuous(limits=c(-1,1.7),breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.3,1,3,10,30))+
  scale_y_continuous(breaks=2:5,limits=c(1.5,5)) +
  scale_color_manual(values=c(color_pd_ecog,color_et_ecog))
ggsave('fig/B02_06_exponent-vs-log-knee_dx.pdf',width=4,height=4)

#------

DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(log_knee>0) %>%
  with(summary(lm(exponent~log_knee)))

DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(log_knee>0) %>%
  with(cor.test(exponent,log_knee,method = 'pearson'))

DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  filter(log_knee>0) %>%
  with(cor.test(exponent,log_knee,method = 'spearman'))

# Spearman's rank correlation rho
# 
# data:  exponent and log_knee
# S = 7.6863e+10, p-value < 2.2e-16
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5627241 

#===========
f_min = 1

aper_curves <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  mutate(dx=factor(dx,levels=c('PD','ET'),ordered=TRUE)) %>%
  mutate(id=paste(subject,session_id,electrode,sep='_')) %>%
  group_by(id,subject,electrode,session_id,electrode_type,dx,dbs_target) %>%
  group_modify(~{ 
    tibble(f=10^seq(0,2.4,0.1)) %>%
      mutate(P= .x$offset * ((10^(.x$log_knee*.x$exponent) + f_min^.x$exponent)/(10^(.x$log_knee*.x$exponent) + f^.x$exponent)))
  })

aper_curves %>%
  ggplot() +
  aes(x=f,y=P,color=dx,group=id)+
  geom_line(alpha=0.03,size=0.5)+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_manual(values=c(color_pd_ecog,color_et_ecog)) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position=c(0.2,0.2), legend.background = element_blank()) +
  coord_cartesian(ylim=c(0.0001,10),xlim=c(1,200),expand=FALSE) 
ggsave('fig/B02_07_aper_curve_overlay.pdf',width=4,height=4)

#===========

DBS_aper_coord_dx

coverage <- DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
  group_by(subject, HCPMMP1_label_1) %>%
  tally(name='n_e') %>%
  group_by(HCPMMP1_label_1) %>%
  summarise(n_e=sum(n_e),n_s=n()) %>%
  arrange(desc(n_s),desc(n_e))

roi <- coverage %>%
  filter(n_s >= 10)

aper_ecog_roi <-  DBS_aper_coord_dx %>%
  filter(electrode_type=='ecog') %>%
  filter(dx%in%c('PD','ET')) %>%
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
ggsave('fig/B02_08_aper_log-knee_vs_anatomy.pdf',width=5,height=4)


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
ggsave('fig/B02_08b_aper_log-knee_vs_anatomy.pdf',width=5,height=4)

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
  select(region1,region,n_s,everything()) %>%
write_tsv('data/aper_param_summary_by_MMP1.tsv')
  


me_mod <- lmer(log_knee ~ region + (1 | subject), 
               data=aper_ecog_roi_agg %>% within(region<-relevel(region,'6v')))
summary(me_mod)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log_knee ~ region + (1 | subject)
#    Data: aper_ecog_roi_agg %>% within(region <- relevel(region, "6v"))
# 
# REML criterion at convergence: -868.4
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.8316 -0.5823  0.0501  0.5262  3.8739 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  subject  (Intercept) 0.005647 0.07514 
#  Residual             0.005341 0.07308 
# Number of obs: 453, groups:  subject, 47
# 
# Fixed effects:
#                Estimate Std. Error         df t value Pr(>|t|)    
# (Intercept)   1.254e+00  1.590e-02  1.450e+02  78.846   <2e-16 ***
# regionp9-46v -1.094e-02  2.476e-02  3.908e+02  -0.442   0.6588    
# region8C     -1.099e-02  1.887e-02  3.885e+02  -0.582   0.5607    
# region8Av    -3.310e-02  1.966e-02  3.885e+02  -1.683   0.0931 .  
# regionIFSp   -1.199e-02  2.394e-02  3.899e+02  -0.501   0.6167    
# region44     -1.158e-02  2.118e-02  3.897e+02  -0.547   0.5847    
# region6r     -4.179e-02  2.117e-02  3.895e+02  -1.974   0.0490 *  
# region55b     1.178e-02  1.839e-02  3.883e+02   0.640   0.5224    
# region4      -2.570e-03  2.153e-02  3.888e+02  -0.119   0.9050    
# region3a     -1.773e-02  2.317e-02  3.892e+02  -0.765   0.4448    
# region3b      3.161e-02  1.918e-02  3.900e+02   1.648   0.1002    
# region1       3.968e-02  1.598e-02  3.869e+02   2.483   0.0135 *  
# region43     -1.281e-04  2.068e-02  3.886e+02  -0.006   0.9951    
# regionOP4    -2.881e-02  1.726e-02  3.872e+02  -1.669   0.0959 .  
# regionPFop   -1.825e-02  1.822e-02  3.888e+02  -1.001   0.3174    
# regionPF      1.092e-03  2.135e-02  3.929e+02   0.051   0.9592    
# regionPSL    -4.418e-02  2.798e-02  3.927e+02  -1.579   0.1152    
# regionSTV    -9.492e-03  2.579e-02  3.929e+02  -0.368   0.7130    
# regionA5     -3.586e-02  1.917e-02  3.892e+02  -1.871   0.0622 .  
# regionA4     -1.124e-02  1.681e-02  3.865e+02  -0.668   0.5042    
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
# 8C - p9-46v == 0   -0.0028331  0.0253873  -0.112   1.0000   
# 8Av - p9-46v == 0  -0.0248175  0.0259427  -0.957   1.0000   
# IFSp - p9-46v == 0 -0.0011821  0.0289391  -0.041   1.0000   
# 44 - p9-46v == 0    0.0059330  0.0269408   0.220   1.0000   
# 6r - p9-46v == 0   -0.0278424  0.0269408  -1.033   1.0000   
# 6v - p9-46v == 0    0.0076282  0.0237264   0.322   1.0000   
# 55b - p9-46v == 0   0.0185986  0.0250805   0.742   1.0000   
# 4 - p9-46v == 0     0.0005177  0.0272559   0.019   1.0000   
# 3a - p9-46v == 0   -0.0122337  0.0284387  -0.430   1.0000   
# 3b - p9-46v == 0    0.0342334  0.0255583   1.339   0.9986   
# 1 - p9-46v == 0     0.0500594  0.0235426   2.126   0.8315   
# 43 - p9-46v == 0    0.0022153  0.0266557   0.083   1.0000   
# OP4 - p9-46v == 0  -0.0190684  0.0243689  -0.782   1.0000   
# PFop - p9-46v == 0 -0.0091862  0.0249424  -0.368   1.0000   
# PF - p9-46v == 0    0.0046271  0.0269408   0.172   1.0000   
# PSL - p9-46v == 0  -0.0332801  0.0318768  -1.044   1.0000   
# STV - p9-46v == 0  -0.0098706  0.0301755  -0.327   1.0000   
# A5 - p9-46v == 0   -0.0269461  0.0255583  -1.054   1.0000   
# A4 - p9-46v == 0   -0.0019613  0.0240966  -0.081   1.0000   
# 8Av - 8C == 0      -0.0219843  0.0211322  -1.040   1.0000   
# IFSp - 8C == 0      0.0016510  0.0247188   0.067   1.0000   
# 44 - 8C == 0        0.0087661  0.0223463   0.392   1.0000   
# 6r - 8C == 0       -0.0250093  0.0223463  -1.119   0.9999   
# 6v - 8C == 0        0.0104614  0.0183437   0.570   1.0000   
# 55b - 8C == 0       0.0214317  0.0200644   1.068   0.9999   
# 4 - 8C == 0         0.0033508  0.0227252   0.147   1.0000   
# 3a - 8C == 0       -0.0094006  0.0241310  -0.390   1.0000   
# 3b - 8C == 0        0.0370665  0.0206585   1.794   0.9595   
# 1 - 8C == 0         0.0528925  0.0181053   2.921   0.2700   
# 43 - 8C == 0        0.0050485  0.0220017   0.229   1.0000   
# OP4 - 8C == 0      -0.0162353  0.0191674  -0.847   1.0000   
# PFop - 8C == 0     -0.0063531  0.0198914  -0.319   1.0000   
# PF - 8C == 0        0.0074603  0.0223463   0.334   1.0000   
# PSL - 8C == 0      -0.0304470  0.0281012  -1.083   0.9999   
# STV - 8C == 0      -0.0070375  0.0261554  -0.269   1.0000   
# A5 - 8C == 0       -0.0241130  0.0206585  -1.167   0.9998   
# A4 - 8C == 0        0.0008718  0.0188201   0.046   1.0000   
# IFSp - 8Av == 0     0.0236353  0.0252888   0.935   1.0000   
# 44 - 8Av == 0       0.0307504  0.0229752   1.338   0.9986   
# 6r - 8Av == 0      -0.0030250  0.0229752  -0.132   1.0000   
# 6v - 8Av == 0       0.0324457  0.0191049   1.698   0.9766   
# 55b - 8Av == 0      0.0434160  0.0207626   2.091   0.8516   
# 4 - 8Av == 0        0.0253351  0.0233439   1.085   0.9999   
# 3a - 8Av == 0       0.0125838  0.0247146   0.509   1.0000   
# 3b - 8Av == 0       0.0590509  0.0213373   2.767   0.3668   
# 1 - 8Av == 0        0.0748768  0.0188761   3.967   0.0114 * 
#   43 - 8Av == 0       0.0270328  0.0226402   1.194   0.9997   
# OP4 - 8Av == 0      0.0057491  0.0198971   0.289   1.0000   
# PFop - 8Av == 0     0.0156312  0.0205955   0.759   1.0000   
# PF - 8Av == 0       0.0294446  0.0229752   1.282   0.9992   
# PSL - 8Av == 0     -0.0084627  0.0286039  -0.296   1.0000   
# STV - 8Av == 0      0.0149468  0.0266948   0.560   1.0000   
# A5 - 8Av == 0      -0.0021286  0.0213373  -0.100   1.0000   
# A4 - 8Av == 0       0.0228562  0.0195627   1.168   0.9998   
# 44 - IFSp == 0      0.0071151  0.0263118   0.270   1.0000   
# 6r - IFSp == 0     -0.0266603  0.0263118  -1.013   1.0000   
# 6v - IFSp == 0      0.0088104  0.0230097   0.383   1.0000   
# 55b - IFSp == 0     0.0197807  0.0244036   0.811   1.0000   
# 4 - IFSp == 0       0.0016998  0.0266343   0.064   1.0000   
# 3a - IFSp == 0     -0.0110515  0.0278435  -0.397   1.0000   
# 3b - IFSp == 0      0.0354156  0.0248943   1.423   0.9970   
# 1 - IFSp == 0       0.0512415  0.0228200   2.245   0.7573   
# 43 - IFSp == 0      0.0033975  0.0260198   0.131   1.0000   
# OP4 - IFSp == 0    -0.0178863  0.0236716  -0.756   1.0000   
# PFop - IFSp == 0   -0.0080041  0.0242615  -0.330   1.0000   
# PF - IFSp == 0      0.0058093  0.0263118   0.221   1.0000   
# PSL - IFSp == 0    -0.0320980  0.0313470  -1.024   1.0000   
# STV - IFSp == 0    -0.0086885  0.0296152  -0.293   1.0000   
# A5 - IFSp == 0     -0.0257640  0.0248943  -1.035   1.0000   
# A4 - IFSp == 0     -0.0007792  0.0233912  -0.033   1.0000   
# 6r - 44 == 0       -0.0337754  0.0240966  -1.402   0.9975   
# 6v - 44 == 0        0.0016953  0.0204397   0.083   1.0000   
# 55b - 44 == 0       0.0126656  0.0219971   0.576   1.0000   
# 4 - 44 == 0        -0.0054153  0.0244484  -0.221   1.0000   
# 3a - 44 == 0       -0.0181666  0.0257604  -0.705   1.0000   
# 3b - 44 == 0        0.0283005  0.0225403   1.256   0.9994   
# 1 - 44 == 0         0.0441264  0.0202260   2.182   0.7980   
# 43 - 44 == 0       -0.0037176  0.0237774  -0.156   1.0000   
# OP4 - 44 == 0      -0.0250014  0.0211821  -1.180   0.9998   
# PFop - 44 == 0     -0.0151192  0.0218394  -0.692   1.0000   
# PF - 44 == 0       -0.0013058  0.0240966  -0.054   1.0000   
# PSL - 44 == 0      -0.0392131  0.0295122  -1.329   0.9988   
# STV - 44 == 0      -0.0158036  0.0276658  -0.571   1.0000   
# A5 - 44 == 0       -0.0328791  0.0225403  -1.459   0.9959   
# A4 - 44 == 0       -0.0078943  0.0208683  -0.378   1.0000   
# 6v - 6r == 0        0.0354707  0.0204397   1.735   0.9708   
# 55b - 6r == 0       0.0464410  0.0219971   2.111   0.8399   
# 4 - 6r == 0         0.0283601  0.0244484   1.160   0.9998   
# 3a - 6r == 0        0.0156088  0.0257604   0.606   1.0000   
# 3b - 6r == 0        0.0620759  0.0225403   2.754   0.3778   
# 1 - 6r == 0         0.0779018  0.0202260   3.852   0.0174 * 
#   43 - 6r == 0        0.0300578  0.0237774   1.264   0.9994   
# OP4 - 6r == 0       0.0087740  0.0211821   0.414   1.0000   
# PFop - 6r == 0      0.0186562  0.0218394   0.854   1.0000   
# PF - 6r == 0        0.0324696  0.0240966   1.347   0.9985   
# PSL - 6r == 0      -0.0054377  0.0295122  -0.184   1.0000   
# STV - 6r == 0       0.0179718  0.0276658   0.650   1.0000   
# A5 - 6r == 0        0.0008964  0.0225403   0.040   1.0000   
# A4 - 6r == 0        0.0258812  0.0208683   1.240   0.9995   
# 55b - 6v == 0       0.0109703  0.0179167   0.612   1.0000   
# 4 - 6v == 0        -0.0071106  0.0208533  -0.341   1.0000   
# 3a - 6v == 0       -0.0198619  0.0223770  -0.888   1.0000   
# 3b - 6v == 0        0.0266052  0.0185796   1.432   0.9967   
# 1 - 6v == 0         0.0424311  0.0156916   2.704   0.4130   
# 43 - 6v == 0       -0.0054129  0.0200625  -0.270   1.0000   
# OP4 - 6v == 0      -0.0266966  0.0169061  -1.579   0.9894   
# PFop - 6v == 0     -0.0168145  0.0177227  -0.949   1.0000   
# PF - 6v == 0       -0.0030011  0.0204397  -0.147   1.0000   
# PSL - 6v == 0      -0.0409084  0.0266102  -1.537   0.9922   
# STV - 6v == 0      -0.0174989  0.0245466  -0.713   1.0000   
# A5 - 6v == 0       -0.0345743  0.0185796  -1.861   0.9428   
# A4 - 6v == 0       -0.0095895  0.0165112  -0.581   1.0000   
# 4 - 55b == 0       -0.0180809  0.0223819  -0.808   1.0000   
# 3a - 55b == 0      -0.0308322  0.0238080  -1.295   0.9991   
# 3b - 55b == 0       0.0156349  0.0202803   0.771   1.0000   
# 1 - 55b == 0        0.0314608  0.0176725   1.780   0.9623   
# 43 - 55b == 0      -0.0163832  0.0216470  -0.757   1.0000   
# OP4 - 55b == 0     -0.0376670  0.0187592  -2.008   0.8911   
# PFop - 55b == 0    -0.0277848  0.0194983  -1.425   0.9969   
# PF - 55b == 0      -0.0139714  0.0219971  -0.635   1.0000   
# PSL - 55b == 0     -0.0518787  0.0278244  -1.865   0.9411   
# STV - 55b == 0     -0.0284692  0.0258577  -1.101   0.9999   
# A5 - 55b == 0      -0.0455447  0.0202803  -2.246   0.7581   
# A4 - 55b == 0      -0.0205599  0.0184041  -1.117   0.9999   
# 3a - 4 == 0        -0.0127514  0.0260897  -0.489   1.0000   
# 3b - 4 == 0         0.0337157  0.0229160   1.471   0.9954   
# 1 - 4 == 0          0.0495417  0.0206439   2.400   0.6465   
# 43 - 4 == 0         0.0016977  0.0241339   0.070   1.0000   
# OP4 - 4 == 0       -0.0195861  0.0215815  -0.908   1.0000   
# PFop - 4 == 0      -0.0097039  0.0222270  -0.437   1.0000   
# PF - 4 == 0         0.0041095  0.0244484   0.168   1.0000   
# PSL - 4 == 0       -0.0337978  0.0298001  -1.134   0.9999   
# STV - 4 == 0       -0.0103883  0.0279728  -0.371   1.0000   
# A5 - 4 == 0        -0.0274638  0.0229160  -1.198   0.9997   
# A4 - 4 == 0        -0.0024790  0.0212735  -0.117   1.0000   
# 3b - 3a == 0        0.0464671  0.0243108   1.911   0.9273   
# 1 - 3a == 0         0.0622930  0.0221820   2.808   0.3414   
# 43 - 3a == 0        0.0144490  0.0254620   0.567   1.0000   
# OP4 - 3a == 0      -0.0068347  0.0230571  -0.296   1.0000   
# PFop - 3a == 0      0.0030474  0.0236624   0.129   1.0000   
# PF - 3a == 0        0.0168608  0.0257604   0.655   1.0000   
# PSL - 3a == 0      -0.0210465  0.0308856  -0.681   1.0000   
# STV - 3a == 0       0.0023630  0.0291264   0.081   1.0000   
# A5 - 3a == 0       -0.0147124  0.0243108  -0.605   1.0000   
# A4 - 3a == 0        0.0102724  0.0227692   0.451   1.0000   
# 1 - 3b == 0         0.0158259  0.0183442   0.863   1.0000   
# 43 - 3b == 0       -0.0320181  0.0221988  -1.442   0.9964   
# OP4 - 3b == 0      -0.0533018  0.0193933  -2.748   0.3819   
# PFop - 3b == 0     -0.0434197  0.0201092  -2.159   0.8129   
# PF - 3b == 0       -0.0296063  0.0225403  -1.313   0.9989   
# PSL - 3b == 0      -0.0675136  0.0282558  -2.389   0.6550   
# STV - 3b == 0      -0.0441041  0.0263214  -1.676   0.9799   
# A5 - 3b == 0       -0.0611795  0.0208683  -2.932   0.2649   
# A4 - 3b == 0       -0.0361947  0.0190500  -1.900   0.9313   
# 43 - 1 == 0        -0.0478440  0.0198447  -2.411   0.6394   
# OP4 - 1 == 0       -0.0691278  0.0166471  -4.153    <0.01 **
#   PFop - 1 == 0      -0.0592456  0.0174758  -3.390   0.0809 . 
# PF - 1 == 0        -0.0454322  0.0202260  -2.246   0.7568   
# PSL - 1 == 0       -0.0833395  0.0264465  -3.151   0.1555   
# STV - 1 == 0       -0.0599300  0.0243689  -2.459   0.6013   
# A5 - 1 == 0        -0.0770054  0.0183442  -4.198    <0.01 **
#   A4 - 1 == 0        -0.0520206  0.0162459  -3.202   0.1357   
# OP4 - 43 == 0      -0.0212837  0.0208183  -1.022   1.0000   
# PFop - 43 == 0     -0.0114016  0.0214867  -0.531   1.0000   
# PF - 43 == 0        0.0024118  0.0237774   0.101   1.0000   
# PSL - 43 == 0      -0.0354955  0.0292522  -1.213   0.9996   
# STV - 43 == 0      -0.0120860  0.0273883  -0.441   1.0000   
# A5 - 43 == 0       -0.0291614  0.0221988  -1.314   0.9989   
# A4 - 43 == 0       -0.0041766  0.0204989  -0.204   1.0000   
# PFop - OP4 == 0     0.0098822  0.0185740   0.532   1.0000   
# PF - OP4 == 0       0.0236955  0.0211821   1.119   0.9999   
# PSL - OP4 == 0     -0.0142117  0.0271846  -0.523   1.0000   
# STV - OP4 == 0      0.0091978  0.0251681   0.365   1.0000   
# A5 - OP4 == 0      -0.0078777  0.0193933  -0.406   1.0000   
# A4 - OP4 == 0       0.0171071  0.0174218   0.982   1.0000   
# PF - PFop == 0      0.0138134  0.0218394   0.632   1.0000   
# PSL - PFop == 0    -0.0240939  0.0276999  -0.870   1.0000   
# STV - PFop == 0    -0.0006844  0.0257237  -0.027   1.0000   
# A5 - PFop == 0     -0.0177599  0.0201092  -0.883   1.0000   
# A4 - PFop == 0      0.0072250  0.0182153   0.397   1.0000   
# PSL - PF == 0      -0.0379073  0.0295122  -1.284   0.9992   
# STV - PF == 0      -0.0144978  0.0276658  -0.524   1.0000   
# A5 - PF == 0       -0.0315732  0.0225403  -1.401   0.9975   
# A4 - PF == 0       -0.0065884  0.0208683  -0.316   1.0000   
# STV - PSL == 0      0.0234095  0.0324919   0.720   1.0000   
# A5 - PSL == 0       0.0063340  0.0282558   0.224   1.0000   
# A4 - PSL == 0       0.0313188  0.0269408   1.163   0.9998   
# A5 - STV == 0      -0.0170755  0.0263214  -0.649   1.0000   
# A4 - STV == 0       0.0079094  0.0249045   0.318   1.0000   
# A4 - A5 == 0        0.0249848  0.0190500   1.312   0.9989   
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
ggsave('fig/B02_08c_aper_exponent_vs_anatomy.pdf',width=5,height=4)

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
ggsave('fig/B02_08c_aper_offset_vs_anatomy.pdf',width=5,height=4)


