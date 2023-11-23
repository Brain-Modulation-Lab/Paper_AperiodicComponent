library(coin)
library(multcomp)
library(lmerTest)
library(tidyverse)
library(ggExtra)

theme_set(theme_bw())

setwd('~/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures')

subjects<-read_tsv('data/subjects.tsv')

DBS_aper <- read_csv('data/DBS_aperiodic_components_all_nonqc.csv') %>%
  select(-X1) %>%
  mutate(session_id=ceiling(session_id)) %>%
  mutate(electrode=tolower(electrode)) %>%
  left_join(subjects %>% select(subject,dx,dbs_target))

unique(DBS_aper$subject)

DBS_aper_coord_dx <- read_tsv('data/A01_DBS_aper_coord_dx.tsv') 

color_et_ecog = '#C4604F' #ET ECoG
color_pd_ecog = '#6F67A6' #PD ECoG

#=============
DBS_aper_agg <- DBS_aper %>%
  filter(electrode_type=='ecog') %>%
  filter(dx=='PD') %>%
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
ggsave('fig/B02_01_exponent-vs-UPDRS.pdf',width=4,height=3)

DBS_aper_agg %>%
  ggplot() +
  aes(x=updrs_score,y=log_knee,color=UPDRS)+
  geom_point() 
  
DBS_aper_agg %>%
  ggplot() +
  aes(x=updrs_score,y=offset,color=UPDRS)+
  geom_point() 

#=============

palette2 <- c('#36A5D1','#C4604F')

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
  scale_y_continuous(limits=c(2.5,4))+
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
# (Intercept)  3.34576    0.09550  35.035   <2e-16 ***
#   dxPD        -0.08029    0.13505  -0.595    0.559    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3167 on 20 degrees of freedom
# Multiple R-squared:  0.01737,	Adjusted R-squared:  -0.03176 
# F-statistic: 0.3535 on 1 and 20 DF,  p-value: 0.5588

oneway_test(exponent ~ factor(dx), data=DBS_aper_agg2, distribution=approximate(nresample=999))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by factor(dx) (ET, PD)
# Z = 0.60391, p-value = 0.5475
# alternative hypothesis: true mu is not equal to 0

summary(lm(log_knee ~ dx, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.22133    0.04443  27.491   <2e-16 ***
#   dxPD        -0.01028    0.06283  -0.164    0.872    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1473 on 20 degrees of freedom
# Multiple R-squared:  0.001337,	Adjusted R-squared:  -0.0486 
# F-statistic: 0.02678 on 1 and 20 DF,  p-value: 0.8716

oneway_test(log_knee ~ factor(dx), data=DBS_aper_agg2, distribution=approximate(nresample=999))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by factor(dx) (ET, PD)
# Z = 0.16758, p-value = 0.8739
# alternative hypothesis: true mu is not equal to 0

summary(lm(offset ~ dx, data=DBS_aper_agg2))
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.726827   0.157571  10.959 6.64e-10 ***
#   dxPD        -0.003964   0.222839  -0.018    0.986    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5226 on 20 degrees of freedom
# Multiple R-squared:  1.582e-05,	Adjusted R-squared:  -0.04998 
# F-statistic: 0.0003164 on 1 and 20 DF,  p-value: 0.986

oneway_test(offset ~ factor(dx), data=DBS_aper_agg2, distribution=approximate(nresample=999))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  offset by factor(dx) (ET, PD)
# Z = 0.018227, p-value = 0.992
# alternative hypothesis: true mu is not equal to 0

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
  filter(n_s >= 7)

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
  arrange(desc(n_r),desc(n_e)) 
  # within(fs_anatomy[HCPMMP1_label_1=='4']<-'precentral') %>%
  # within(fs_anatomy[HCPMMP1_label_1=='3b']<-'postcentral')
  
aper_ecog_roi <- within(aper_ecog_roi,fs_anatomy[HCPMMP1_label_1=='4']<-'precentral')
aper_ecog_roi <- within(aper_ecog_roi,fs_anatomy[HCPMMP1_label_1=='3b']<-'postcentral')



#ordering regions by y mni coordinate
regions_coords <- aper_ecog_roi %>%
  group_by(HCPMMP1_label_1) %>%
  summarise(x=mean(mni_nonlinear_x),y=mean(mni_nonlinear_y),z=mean(mni_nonlinear_z),
            fs_anatomy=names(sort(table(fs_anatomy), decreasing = TRUE)[1])) %>%
  group_by(fs_anatomy) %>%
  mutate(fs_anatomy_y=mean(y)) %>%
  ungroup %>%
  arrange(fs_anatomy_y,y)

aper_ecog_roi_agg <- aper_ecog_roi %>%
  filter(log_knee > 1) %>%
  group_by(subject,HCPMMP1_label_1) %>%
  summarise(exponent = median(exponent), log_knee = median(log_knee), offset = median(offset),
            x=mean(mni_nonlinear_x),y=mean(mni_nonlinear_y),z=mean(mni_nonlinear_z)) %>%
  group_by(subject) %>% 
  mutate(exponent_s = exponent - median(exponent), log_knee_s = log_knee - median(log_knee), offset_s = offset - median(offset)) %>%
  ungroup() %>%
  mutate(region = factor(HCPMMP1_label_1,levels=regions_coords$HCPMMP1_label_1)) %>%
  left_join(regions_coords %>% dplyr::select(HCPMMP1_label_1,fs_anatomy)) 

aper_ecog_roi_agg %>%
  ggplot() + 
  aes(x=region,y=exponent,color=subject,group=subject)+
  geom_point() + 
  geom_line() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1))

aper_ecog_roi_agg %>%
  ggplot() + 
  aes(x=region,y=exponent_s,color=fs_anatomy)+
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
  mutate(region=reorder(region, desc(region))) %>%
  ggplot() + 
  aes(x=region,y=log_knee_s,color=fs_anatomy)+
  geom_boxplot() + 
  scale_color_discrete(breaks=rev(unique(regions_coords$fs_anatomy))) +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1)) +
  theme(legend.position='top',legend.title = element_blank(),
        axis.title.x=element_blank(),legend.text=element_text(size=6)) +
  labs(y="Log10(Knee frequency)\n region effect")
ggsave('fig/B02_08_aper_log-knee_vs_anatomy.pdf',width=5,height=4)

me_mod <- lmer(log_knee ~ region + (1 | subject), 
               data=aper_ecog_roi_agg %>% within(region<-relevel(region,5)))
summary(me_mod)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# subject  (Intercept) 0.007778 0.08819 
# Residual             0.002761 0.05254 
# Number of obs: 141, groups:  subject, 20
# 
# Fixed effects:
#   Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)   1.20866    0.02624  47.99490  46.059  < 2e-16 ***
#   regionPSL     0.03425    0.03051 109.11309   1.123 0.264071    
# regionPF      0.03379    0.02783 109.99550   1.214 0.227288    
# regionSTV     0.07028    0.02918 109.56937   2.409 0.017673 *  
#   regionA4      0.02814    0.02219 108.03115   1.268 0.207486    
# regionPFop    0.00479    0.02443 108.55237   0.196 0.844944    
# regionOP4     0.01425    0.02344 108.14641   0.608 0.544585    
# region1       0.07538    0.02108 108.33122   3.576 0.000523 ***
#   region4       0.09508    0.02783 108.09504   3.417 0.000895 ***
#   region3b      0.06501    0.02542 109.22315   2.558 0.011901 *  
#   region55b     0.05000    0.02378 108.51431   2.103 0.037796 *  
#   region43      0.05529    0.02563 108.35046   2.157 0.033196 *  
#   region6v      0.04502    0.02149 108.20145   2.095 0.038526 *  
#   region8Av     0.05457    0.02586 108.58416   2.110 0.037117 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

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
# PF - PSL == 0    0.0032946  0.0298515   0.110   1.0000   
# STV - PSL == 0   0.0370069  0.0308706   1.199   0.9948   
# A4 - PSL == 0   -0.0002356  0.0265606  -0.009   1.0000   
# A5 - PSL == 0   -0.0335401  0.0279235  -1.201   0.9948   
# PFop - PSL == 0 -0.0230922  0.0279235  -0.827   0.9999   
# OP4 - PSL == 0  -0.0109734  0.0274972  -0.399   1.0000   
# 1 - PSL == 0     0.0464268  0.0256244   1.812   0.8608   
# 4 - PSL == 0     0.0586248  0.0308706   1.899   0.8161   
# 3b - PSL == 0    0.0321894  0.0284359   1.132   0.9970   
# 55b - PSL == 0   0.0141734  0.0274972   0.515   1.0000   
# 43 - PSL == 0    0.0304819  0.0290637   1.049   0.9986   
# 6v - PSL == 0    0.0194272  0.0259365   0.749   1.0000   
# 8Av - PSL == 0   0.0250893  0.0290637   0.863   0.9998   
# STV - PF == 0    0.0337123  0.0283633   1.189   0.9952   
# A4 - PF == 0    -0.0035302  0.0235997  -0.150   1.0000   
# A5 - PF == 0    -0.0368347  0.0251238  -1.466   0.9695   
# PFop - PF == 0  -0.0263867  0.0251238  -1.050   0.9986   
# OP4 - PF == 0   -0.0142680  0.0246491  -0.579   1.0000   
# 1 - PF == 0      0.0431322  0.0225408   1.914   0.8080   
# 4 - PF == 0      0.0553302  0.0283633   1.951   0.7864   
# 3b - PF == 0     0.0288948  0.0256921   1.125   0.9972   
# 55b - PF == 0    0.0108788  0.0246491   0.441   1.0000   
# 43 - PF == 0     0.0271873  0.0263852   1.030   0.9988   
# 6v - PF == 0     0.0161326  0.0228950   0.705   1.0000   
# 8Av - PF == 0    0.0217947  0.0263852   0.826   0.9999   
# A4 - STV == 0   -0.0372425  0.0248762  -1.497   0.9638   
# A5 - STV == 0   -0.0705469  0.0263265  -2.680   0.2931   
# PFop - STV == 0 -0.0600990  0.0263265  -2.283   0.5606   
# OP4 - STV == 0  -0.0479803  0.0258739  -1.854   0.8404   
# 1 - STV == 0     0.0094199  0.0238741   0.395   1.0000   
# 4 - STV == 0     0.0216179  0.0294340   0.734   1.0000   
# 3b - STV == 0   -0.0048175  0.0268694  -0.179   1.0000   
# 55b - STV == 0  -0.0228335  0.0258739  -0.882   0.9998   
# 43 - STV == 0   -0.0065250  0.0275329  -0.237   1.0000   
# 6v - STV == 0   -0.0175797  0.0242088  -0.726   1.0000   
# 8Av - STV == 0  -0.0119176  0.0275329  -0.433   1.0000   
# A5 - A4 == 0    -0.0333045  0.0211082  -1.578   0.9461   
# PFop - A4 == 0  -0.0228566  0.0211082  -1.083   0.9981   
# OP4 - A4 == 0   -0.0107378  0.0205409  -0.523   1.0000   
# 1 - A4 == 0      0.0466624  0.0179567   2.599   0.3419   
# 4 - A4 == 0      0.0588604  0.0248762   2.366   0.4996   
# 3b - A4 == 0     0.0324250  0.0217815   1.489   0.9654   
# 55b - A4 == 0    0.0144090  0.0205409   0.701   1.0000   
# 43 - A4 == 0     0.0307174  0.0225950   1.359   0.9837   
# 6v - A4 == 0     0.0196628  0.0183993   1.069   0.9983   
# 8Av - A4 == 0    0.0253249  0.0225950   1.121   0.9973   
# PFop - A5 == 0   0.0104479  0.0227994   0.458   1.0000   
# OP4 - A5 == 0    0.0225667  0.0222753   1.013   0.9990   
# 1 - A5 == 0      0.0799668  0.0199173   4.015    <0.01 **
#   4 - A5 == 0      0.0921649  0.0263265   3.501   0.0374 * 
#   3b - A5 == 0     0.0657294  0.0234242   2.806   0.2277   
# 55b - A5 == 0    0.0477135  0.0222753   2.142   0.6615   
# 43 - A5 == 0     0.0640219  0.0241825   2.647   0.3133   
# 6v - A5 == 0     0.0529673  0.0203173   2.607   0.3381   
# 8Av - A5 == 0    0.0586293  0.0241825   2.424   0.4581   
# OP4 - PFop == 0  0.0121188  0.0222753   0.544   1.0000   
# 1 - PFop == 0    0.0695189  0.0199173   3.490   0.0389 * 
#   4 - PFop == 0    0.0817170  0.0263265   3.104   0.1140   
# 3b - PFop == 0   0.0552815  0.0234242   2.360   0.5054   
# 55b - PFop == 0  0.0372656  0.0222753   1.673   0.9173   
# 43 - PFop == 0   0.0535740  0.0241825   2.215   0.6087   
# 6v - PFop == 0   0.0425194  0.0203173   2.093   0.6957   
# 8Av - PFop == 0  0.0481814  0.0241825   1.992   0.7607   
# 1 - OP4 == 0     0.0574002  0.0193151   2.972   0.1584   
# 4 - OP4 == 0     0.0695982  0.0258739   2.690   0.2896   
# 3b - OP4 == 0    0.0431628  0.0229143   1.884   0.8248   
# 55b - OP4 == 0   0.0251468  0.0217384   1.157   0.9963   
# 43 - OP4 == 0    0.0414553  0.0236889   1.750   0.8885   
# 6v - OP4 == 0    0.0304006  0.0197273   1.541   0.9548   
# 8Av - OP4 == 0   0.0360627  0.0236889   1.522   0.9589   
# 4 - 1 == 0       0.0121980  0.0238741   0.511   1.0000   
# 3b - 1 == 0     -0.0142374  0.0206296  -0.690   1.0000   
# 55b - 1 == 0    -0.0322534  0.0193151  -1.670   0.9187   
# 43 - 1 == 0     -0.0159449  0.0214867  -0.742   1.0000   
# 6v - 1 == 0     -0.0269996  0.0170200  -1.586   0.9438   
# 8Av - 1 == 0    -0.0213375  0.0214867  -0.993   0.9992   
# 3b - 4 == 0     -0.0264354  0.0268694  -0.984   0.9993   
# 55b - 4 == 0    -0.0444514  0.0258739  -1.718   0.9011   
# 43 - 4 == 0     -0.0281429  0.0275329  -1.022   0.9989   
# 6v - 4 == 0     -0.0391976  0.0242088  -1.619   0.9347   
# 8Av - 4 == 0    -0.0335355  0.0275329  -1.218   0.9940   
# 55b - 3b == 0   -0.0180160  0.0229143  -0.786   0.9999   
# 43 - 3b == 0    -0.0017075  0.0247724  -0.069   1.0000   
# 6v - 3b == 0    -0.0127622  0.0210160  -0.607   1.0000   
# 8Av - 3b == 0   -0.0071001  0.0247724  -0.287   1.0000   
# 43 - 55b == 0    0.0163085  0.0236889   0.688   1.0000   
# 6v - 55b == 0    0.0052538  0.0197273   0.266   1.0000   
# 8Av - 55b == 0   0.0109159  0.0236889   0.461   1.0000   
# 6v - 43 == 0    -0.0110546  0.0218580  -0.506   1.0000   
# 8Av - 43 == 0   -0.0053926  0.0254906  -0.212   1.0000   
# 8Av - 6v == 0    0.0056621  0.0218580   0.259   1.0000   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

