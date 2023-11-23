library(tidyverse)
library(GGally)
library(scales)


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



##### Sensitivity analysis all electrodes
range_exponent = c(0.5, 5.5)
range_offset = c(-0.5, 6)
range_log_knee = c(-1, 2)
bw_exponent = 0.05*(range_exponent[2]-range_exponent[1])
bw_offset = 0.05*(range_offset[2]-range_offset[1])
bw_log_knee = 0.05*(range_log_knee[2]-range_log_knee[1])

db_fig_A03_01 <- db %>%
  group_by(subject,run,label) %>%
  mutate(is_ref = sum(freq_range=="1-250hz")==1) %>%
  filter(is_ref) %>%
  mutate(offset_ref = offset[freq_range=="1-250hz"],
         exponent_ref = exponent[freq_range=="1-250hz"],
         log_knee_ref = log_knee[freq_range=="1-250hz"]) %>%
  mutate(freq_range=factor(freq_range,levels=c('1-250hz','1-200hz','1-150hz'))) %>%
  write_tsv('fig/A03_01_senitivity-analysis.tsv')

db_fig_A03_01_percent <- db_fig_A03_01 %>% 
  group_by(freq_range) %>% 
  summarise(f_offset = sum(abs(offset - offset_ref)<bw_offset)/length(offset),
            f_exponent = sum(abs(exponent - exponent_ref)<bw_exponent)/length(exponent),
            f_log_knee = sum(abs(log_knee - log_knee_ref)<bw_log_knee)/length(log_knee)) %>%
  mutate(label_offset=str_c(signif(f_offset*100,3),'%'),
         label_exponent=str_c(signif(f_exponent*100,3),'%'),
         label_log_knee=str_c(signif(f_log_knee*100,3),'%')) %>%
  write_tsv('fig/A03_01_senitivity-analysis-percent.tsv')

db_fig_A03_01 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = after_stat(count)),binwidth=bw_offset,alpha=0.9,alpha=0.9) +
  #scale_fill_A03_03 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,6))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=6,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_01a_senitivity-analysis-offset-hex-bins.pdf',width=7,height=3)


db_fig_A03_01 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_point(alpha=0.5,color='#08519C',size=0.5) +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,4.5),ylim=c(-0.5,4.5))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=4.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_01a_senitivity-analysis-offset-points.pdf',width=6,height=3)

db_fig_A03_01 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_point(alpha=0.5,color='#08519C',size=0.5) +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(1,4.5),ylim=c(1,4.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=1,y=4,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_01a_senitivity-analysis-exponent-points.pdf',width=6,height=3)


db_fig_A03_01 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_point(alpha=0.5,color='#08519C',size=0.5) +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-1,2),ylim=c(-1,2))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=2,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_01a_senitivity-analysis-log_knee-points.pdf',width=6,height=3)



#color by type

color_ep_ctx = '#C4604F' 
color_ep_thal = '#36A5D1' #ET VIM


db_fig_A03_01 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  geom_point(aes(color=type),alpha=0.5,size=0.5) +
  scale_color_manual(values=c(color_ep_ctx,color_ep_thal)) +
  coord_equal(xlim=c(-0.5,4.5),ylim=c(-0.5,4.5))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=4.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_02_senitivity-analysis-offset-points.pdf',width=6,height=3)


db_fig_A03_01 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  geom_point(aes(color=type),alpha=0.5,size=0.5) +
  scale_color_manual(values=c(color_ep_ctx,color_ep_thal)) +
  coord_equal(xlim=c(1,4.5),ylim=c(1,4.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=1,y=4,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_02_senitivity-analysis-exponent-points.pdf',width=6,height=3)


db_fig_A03_01 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  geom_point(aes(color=type),alpha=0.5,size=0.5) +
  scale_color_manual(values=c(color_ep_ctx,color_ep_thal)) +
  coord_equal(xlim=c(-1,1.5),ylim=c(-1,1.5))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=1.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_01_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_02_senitivity-analysis-log_knee-points.pdf',width=6,height=3)


#####

scale_fill_A03_03 <-
  scale_fill_gradientn(limits=c(1,30),oob=scales::squish,
                       colours=c('#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6','#2171B5','#08519C','#08306B')) 

##### only ctx

db_fig_A03_03 <- db_fig_A03_01 %>% 
  filter(type=='ctx')

db_fig_A03_03_percent <- db_fig_A03_03 %>% 
  group_by(freq_range) %>% 
  summarise(f_offset = sum(abs(offset - offset_ref)<bw_offset)/length(offset),
            f_exponent = sum(abs(exponent - exponent_ref)<bw_exponent)/length(exponent),
            f_log_knee = sum(abs(log_knee - log_knee_ref)<bw_log_knee)/length(log_knee)) %>%
  mutate(label_offset=str_c(signif(f_offset*100,3),'%'),
         label_exponent=str_c(signif(f_exponent*100,3),'%'),
         label_log_knee=str_c(signif(f_log_knee*100,3),'%')) 

db_fig_A03_03 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_offset,alpha=0.9) +
  scale_fill_A03_03 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,6),ylim=c(-0.5,6))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=6,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_03_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_03_senitivity-analysis-offset-ctx.pdf',width=7,height=3)

db_fig_A03_03 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_exponent,alpha=0.9) +
  scale_fill_A03_03 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(0.5,5.5),ylim=c(0.5,5.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=0.5,y=5.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_03_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_03_senitivity-analysis-exponent-ctx.pdf',width=7,height=3)

db_fig_A03_03 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_log_knee,alpha=0.9) +
  scale_fill_A03_03 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-1,2),ylim=c(-1,2))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=2,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_03_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_03_senitivity-analysis-log_knee-ctx.pdf',width=7,height=3)




##### only thal
db_fig_A03_04 <- db_fig_A03_01 %>% 
  filter(type=='thal')

db_fig_A03_04_percent <- db_fig_A03_04 %>% 
  group_by(freq_range) %>% 
  summarise(f_offset = sum(abs(offset - offset_ref)<bw_offset)/length(offset),
            f_exponent = sum(abs(exponent - exponent_ref)<bw_exponent)/length(exponent),
            f_log_knee = sum(abs(log_knee - log_knee_ref)<bw_log_knee)/length(log_knee)) %>%
  mutate(label_offset=str_c(signif(f_offset*100,3),'%'),
         label_exponent=str_c(signif(f_exponent*100,3),'%'),
         label_log_knee=str_c(signif(f_log_knee*100,3),'%')) 

scale_fill_A03_04 <-
  scale_fill_gradientn(limits=c(0,30),oob=scales::squish,
                       colours=c('#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6','#2171B5','#08519C','#08306B')) 

db_fig_A03_04 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_offset,alpha=0.9) +
  scale_fill_A03_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,6),ylim=c(-0.5,6))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=6,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_04_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_04_senitivity-analysis-offset-thal.pdf',width=7,height=3)


db_fig_A03_04 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_exponent,alpha=0.9) +
  scale_fill_A03_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(0.5,5.5),ylim=c(0.5,5.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=0.5,y=5.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_04_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_04_senitivity-analysis-exponent-thal.pdf',width=7,height=3)

db_fig_A03_04 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_log_knee,alpha=0.9) +
  scale_fill_A03_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-1,2),ylim=c(-1,2))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=2,vjust=1,hjust=0,color='gray60',
            data=db_fig_A03_04_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/A03_04_senitivity-analysis-log_knee-thal.pdf',width=7,height=3)



db_fig_A03_04 %>%
  filter(freq_range=='1-150hz') %>%
  mutate(has_knee_ref = log_knee_ref > 0, has_knee = log_knee >0) %>%
  group_by(has_knee_ref,has_knee) %>%
  summarise(ne=n(),ns=length(unique(subject))) 

#   has_knee_ref has_knee    ne    ns
#   <lgl>        <lgl>    <int> <int>
# 1 FALSE        FALSE       41     6
# 2 FALSE        TRUE        61     8
# 3 TRUE         FALSE        4     2
# 4 TRUE         TRUE        30     7

