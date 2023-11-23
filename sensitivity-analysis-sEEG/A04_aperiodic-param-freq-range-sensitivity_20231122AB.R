library(tidyverse)
library(GGally)
library(scales)
library(coin)

theme_set(theme_bw())

PATH_ANALISYS <- '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-11-22-FOOOF-sensitivity-analysis-sEEG'
setwd(PATH_ANALISYS)

get_type <- function(is_cortical,is_thalamic){
  if(is_cortical)
    return('ctx')
  else if(is_thalamic)
    return('thal')
  else
    return('none')
}

elecloc <- read_csv('data/electrode-localization.csv') %>%
  select(-'...1') %>%
  rename(c(subject=Subject)) %>%
  rowwise() %>%
  mutate(type=get_type(is_cortical,is_thalamic)) %>%
  select(-is_thalamic,-is_cortical)

db <- read_csv('data/spectral_sensitivity_ap-params.csv') %>%
  select(-'...1') %>%
  rename(c(subject=Subject)) %>%
  left_join(elecloc) %>%
  filter(type %in% c('thal','ctx')) %>%
  mutate(offset = offset + 12) #unit conversion from V^2 to uV^2


#######

color_ep_cm = '#36A5D1' #EP CM
color_ep_seeg = '#C4604F' #EP sEEG

se <- function(x) sd(x)/sqrt(length(x))

sEEG_aper_agg2 <- db %>%
  mutate(freq_range=factor(freq_range,levels=c('1-250hz','1-200hz','1-150hz'))) %>%
  mutate(type=factor(type,levels=c('thal','ctx'))) %>%
  group_by(subject,type,freq_range) %>%
  summarise(se_offset=se(offset),se_log_knee=se(log_knee),se_exponent=se(exponent),
            sd_offset=sd(offset),sd_log_knee=sd(log_knee),sd_exponent=sd(exponent),
            offset=median(offset),log_knee=median(log_knee),exponent=median(exponent)) 

sEEG_aper_agg2 %>%
  ggplot() +
  aes(y=exponent,ymin=exponent-sd_exponent,ymax=exponent+sd_exponent,x=type,color=type)+
  #geom_errorbar(position=position_dodge2(width = 0.5)) +
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  #scale_y_continuous(limits=c(2.5,4))+
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  labs(y='Exponent')+
  facet_wrap(~freq_range)
ggsave('fig/A04_01_exponent-vs-electrode-type_EP.pdf',width=6,height=2.5)

sEEG_aper_agg2 %>%
  ggplot() +
  aes(y=log_knee,ymin=log_knee-sd_log_knee,ymax=log_knee+sd_log_knee,x=type,color=type)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  geom_hline(yintercept=0,color='gray',linetype=2) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Knee Frequency (Hz)')+
  facet_wrap(~freq_range)
ggsave('fig/A04_02_log-knee-vs-electrode-type_EP.pdf',width=6,height=2.5)


sEEG_aper_agg2 %>%
  ggplot() +
  aes(y=offset,ymin=offset-sd_offset,ymax=offset+sd_offset,x=type,color=type)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  #scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  #coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Offset')+
  facet_wrap(~freq_range)
ggsave('fig/A04_03_offset-vs-electrode-type_EP.pdf',width=6,height=2.5)

###







sEEG_aper_agg2 %>%
  filter(freq_range=='1-150hz') %>%
  ggplot() +
  aes(y=exponent,ymin=exponent-sd_exponent,ymax=exponent+sd_exponent,x=type,color=type)+
  #geom_errorbar(position=position_dodge2(width = 0.5)) +
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  #scale_y_continuous(limits=c(2.5,4))+
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  labs(y='Exponent')+
  facet_wrap(~freq_range)
ggsave('fig/A05_01_exponent-vs-electrode-type_EP_1-150hz.pdf',width=2.5,height=2.5)

sEEG_aper_agg2 %>%
  filter(freq_range=='1-150hz') %>%
  ggplot() +
  aes(y=log_knee,ymin=log_knee-sd_log_knee,ymax=log_knee+sd_log_knee,x=type,color=type)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  geom_hline(yintercept=0,color='gray',linetype=2) +
  scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Knee Frequency (Hz)')+
  facet_wrap(~freq_range)
ggsave('fig/A05_02_log-knee-vs-electrode-type_EP_1-150hz.pdf',width=2.5,height=2.5)


sEEG_aper_agg2 %>%
  filter(freq_range=='1-150hz') %>%
  ggplot() +
  aes(y=offset,ymin=offset-sd_offset,ymax=offset+sd_offset,x=type,color=type)+
  geom_line(aes(group = subject),position=position_dodge2(width = 0.75),color='gray') +
  geom_pointrange(position=position_dodge2(width = 0.75),size=0.2) +
  #scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5),labels=c(0.1,0.31,1,3.1,10,31))+
  theme(legend.position='none',panel.grid = element_blank(), axis.title.x=element_blank()) +
  scale_color_manual(values=c(color_ep_cm,color_ep_seeg)) +
  #coord_cartesian(ylim=c(-1.1,1.5)) +
  labs(y='Offset')+
  facet_wrap(~freq_range)
ggsave('fig/A05_03_offset-vs-electrode-type_EP_1-150hz.pdf',width=2.5,height=2.5)



###

oneway_test_pvalue<- function(d,y){
  e<-oneway_test(exponent ~ factor(type), distribution=approximate(nresample=9999),data=d)
  o<-oneway_test(offset ~ factor(type), distribution=approximate(nresample=9999),data=d)
  k<-oneway_test(log_knee ~ factor(type), distribution=approximate(nresample=9999),data=d)
  e<-max(as.numeric(pvalue(e)),1e-4)
  o<-max(as.numeric(pvalue(o)),1e-4)
  k<-max(as.numeric(pvalue(k)),1e-4)
  tibble(exponent=e,offset=o,log_knee=k)
}

db %>%
  group_by(freq_range,subject) %>%
  group_modify(oneway_test_pvalue) %>%
  pivot_longer(cols=c(exponent,offset,log_knee),values_to='pvalue',names_to = 'param') %>%
  ungroup() %>%
  mutate(pvalue_fdr = p.adjust(pvalue,method="fdr")) %>% 
  write_tsv('fig/A04_04_pvalues-fdr.tsv') %>%
  filter(pvalue_fdr > 0.05)

db %>%
  filter(freq_range=='1-150hz') %>%
  group_by(subject) %>%
  group_modify(oneway_test_pvalue) %>%
  pivot_longer(cols=c(exponent,offset,log_knee),values_to='pvalue',names_to = 'param') %>%
  ungroup() %>%
  mutate(pvalue_fdr = p.adjust(pvalue,method="fdr")) %>% 
  filter(pvalue_fdr > 0.05)



oneway_test(log_knee ~ factor(type), distribution=approximate(nresample=9999),
            data=sEEG_aper_agg2 %>% filter(freq_range=='1-150hz'))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by factor(type) (thal, ctx)
# Z = -2.7888, p-value = 0.0007001
# alternative hypothesis: true mu is not equal to 0

oneway_test(exponent ~ factor(type), distribution=approximate(nresample=9999),
            data=sEEG_aper_agg2 %>% filter(freq_range=='1-150hz'))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by factor(type) (thal, ctx)
# Z = -3.1225, p-value = 3e-04
# alternative hypothesis: true mu is not equal to 0

oneway_test(offset ~ factor(type), distribution=approximate(nresample=9999),
            data=sEEG_aper_agg2 %>% filter(freq_range=='1-150hz'))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  offset by factor(type) (thal, ctx)
# Z = -1.2571, p-value = 0.2093
# alternative hypothesis: true mu is not equal to 0


sEEG_aper_agg2_for_test <- sEEG_aper_agg2 %>%
  filter(freq_range=='1-150hz') %>%
  group_by(subject) %>%
  select(subject, type,offset,log_knee,exponent) %>%
  pivot_wider(names_from=type, values_from=c(offset,log_knee,exponent))


sEEG_aper_agg2_for_test %>%
  with(wilcox.test(offset_thal,offset_ctx, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  offset_thal and offset_ctx
# V = 3, p-value = 0.03906
# alternative hypothesis: true location shift is not equal to 0


sEEG_aper_agg2_for_test %>%
  with(wilcox.test(exponent_thal,exponent_ctx, paired=TRUE))
# Wilcoxon signed rank exact test
# 
# data:  exponent_thal and exponent_ctx
# V = 0, p-value = 0.007813
# alternative hypothesis: true location shift is not equal to 0

sEEG_aper_agg2_for_test %>%
  with(wilcox.test(log_knee_thal,log_knee_ctx, paired=TRUE))

# Wilcoxon signed rank exact test
# 
# data:  log_knee_thal and log_knee_ctx
# V = 0, p-value = 0.007813
# alternative hypothesis: true location shift is not equal to 0



#250Hz

oneway_test(log_knee ~ factor(type), distribution=approximate(nresample=9999),
            data=sEEG_aper_agg2 %>% filter(freq_range=='1-250hz'))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  log_knee by factor(type) (thal, ctx)
# Z = -3.5218, p-value = 4e-04
# alternative hypothesis: true mu is not equal to 0

oneway_test(exponent ~ factor(type), distribution=approximate(nresample=9999),
            data=sEEG_aper_agg2 %>% filter(freq_range=='1-250hz'))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  exponent by factor(type) (thal, ctx)
# Z = -3.4296, p-value = 1e-04
# alternative hypothesis: true mu is not equal to 0

oneway_test(offset ~ factor(type), distribution=approximate(nresample=9999),
            data=sEEG_aper_agg2 %>% filter(freq_range=='1-250hz'))
# Approximative Two-Sample Fisher-Pitman Permutation Test
# 
# data:  offset by factor(type) (thal, ctx)
# Z = -1.3323, p-value = 0.1893
# alternative hypothesis: true mu is not equal to 0