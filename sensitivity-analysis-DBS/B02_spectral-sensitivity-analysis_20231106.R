library(tidyverse)
library(GGally)
library(scales)
theme_set(theme_bw())


PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-10-31-FOOOF-sensitivity-analysis'
setwd(PATH_ANALYSIS)

db <- read_csv('fooof_outputs/spectral_sensitivity_ap-params.csv') %>%
  select(-'...1') %>%
  filter(str_detect(electrode,'ecog|dbs_L')) %>%
  filter(R2 > 0.95) %>%
  mutate(electrode_type = factor(str_split_i(electrode,"_",1),levels=c('ecog','dbs')))

db_agg <- db %>%
  group_by(subject,freq_range,electrode,electrode_type) %>%
  summarise(offset = median(offset), 
            log_knee=median(log_knee),
            exponent=median(exponent),
            R2=median(R2),
            error=median(error))


db %>%
  pivot_wider(id_cols=c(subject,session_id,electrode,electrode_type),names_from=freq_range,values_from=offset) %>%
  ggpairs(columns=5:10,title="offset",axisLabels="show",mapping=aes(alpha=0.01,color=electrode_type)) +
  scale_alpha_identity() +
  geom_abline(slope=1,intercept=0,color="gray")
ggsave('fig/B02_01_ggpairs-offset.png',width=9,height=9)

db %>%
  pivot_wider(id_cols=c(subject,session_id,electrode,electrode_type),names_from=freq_range,values_from=log_knee) %>%
  ggpairs(columns=5:10,title="log_knee",axisLabels="show",mapping=aes(alpha=0.01,color=electrode_type)) +
  scale_alpha_identity() +
  geom_abline(slope=1,intercept=0,color="gray")
ggsave('fig/B02_02_ggpairs-log-knee.png',width=9,height=9)

db %>%
  pivot_wider(id_cols=c(subject,session_id,electrode,electrode_type),names_from=freq_range,values_from=exponent) %>%
  ggpairs(columns=5:10,title="exponent",axisLabels="show",mapping=aes(alpha=0.01,color=electrode_type)) +
  scale_alpha_identity() +
  geom_abline(slope=1,intercept=0,color="gray")
ggsave('fig/B02_03_ggpairs-exponent.png',width=9,height=9)


##### Sensitivity analysis all electrodes
range_exponent = c(0.5, 5.5)
range_offset = c(-0.5, 6)
range_log_knee = c(-1, 2)
bw_exponent = 0.05*(range_exponent[2]-range_exponent[1])
bw_offset = 0.05*(range_offset[2]-range_offset[1])
bw_log_knee = 0.05*(range_log_knee[2]-range_log_knee[1])


db_fig_B02_04 <- db %>%
  group_by(subject,session_id,electrode) %>%
  mutate(is_ref = sum(freq_range=="1-250hz")==1) %>%
  filter(is_ref) %>%
  mutate(offset_ref = offset[freq_range=="1-250hz"],
         exponent_ref = exponent[freq_range=="1-250hz"],
         log_knee_ref = log_knee[freq_range=="1-250hz"]) %>%
  mutate(freq_range=factor(freq_range,levels=c('1-250hz','1-200hz','1-150hz','5-250hz','30-250hz','1-12--30-250hz'))) %>%
  write_tsv('fig/B02_04_senitivity-analysis.tsv')

db_fig_B02_04_percent <- db_fig_B02_04 %>% 
  group_by(freq_range) %>% 
  summarise(f_offset = sum(abs(offset - offset_ref)<bw_offset)/length(offset),
            f_exponent = sum(abs(exponent - exponent_ref)<bw_exponent)/length(exponent),
            f_log_knee = sum(abs(log_knee - log_knee_ref)<bw_log_knee)/length(log_knee)) %>%
  mutate(label_offset=str_c(signif(f_offset*100,3),'%'),
         label_exponent=str_c(signif(f_exponent*100,3),'%'),
         label_log_knee=str_c(signif(f_log_knee*100,3),'%')) %>%
  write_tsv('fig/B02_04_senitivity-analysis-percent.tsv')

scale_fill_B02_04 <-
  scale_fill_gradientn(limits=c(1,1000),oob=scales::squish,
                       colours=c('#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6','#2171B5','#08519C','#08306B')) 

db_fig_B02_04 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_offset,alpha=0.9,alpha=0.9) +
  scale_fill_B02_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,6))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=6,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_04_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_04a_senitivity-analysis-offset.pdf',width=10,height=3)

db_fig_B02_04 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_exponent,alpha=0.9) +
  scale_fill_B02_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(0.5,5.5),ylim=c(0.5,5.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=0.5,y=5.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_04_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_04b_senitivity-analysis-exponent.pdf',width=10,height=3)

db_fig_B02_04 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_log_knee,alpha=0.9) +
  scale_fill_B02_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-1,2),ylim=c(-1,2))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=2,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_04_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_04c_senitivity-analysis-log_knee.pdf',width=10,height=3)


##### only ECoG
db_fig_B02_05 <- db_fig_B02_04 %>% 
  filter(str_starts(electrode,'ecog'))

db_fig_B02_05_percent <- db_fig_B02_05 %>% 
  filter(str_starts(electrode,'ecog')) %>%
  group_by(freq_range) %>% 
  summarise(f_offset = sum(abs(offset - offset_ref)<bw_offset)/length(offset),
            f_exponent = sum(abs(exponent - exponent_ref)<bw_exponent)/length(exponent),
            f_log_knee = sum(abs(log_knee - log_knee_ref)<bw_log_knee)/length(log_knee)) %>%
  mutate(label_offset=str_c(signif(f_offset*100,3),'%'),
         label_exponent=str_c(signif(f_exponent*100,3),'%'),
         label_log_knee=str_c(signif(f_log_knee*100,3),'%')) 

db_fig_B02_05 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_offset,alpha=0.9) +
  scale_fill_B02_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,6),ylim=c(-0.5,6))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=6,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_05_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_05a_senitivity-analysis-offset-ecog.pdf',width=10,height=3)


db_fig_B02_05 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_exponent,alpha=0.9) +
  scale_fill_B02_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(0.5,5.5),ylim=c(0.5,5.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=0.5,y=5.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_05_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_05b_senitivity-analysis-exponent-ecog.pdf',width=10,height=3)

db_fig_B02_05 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_log_knee,alpha=0.9) +
  scale_fill_B02_04 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-1,2),ylim=c(-1,2))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=2,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_05_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_05c_senitivity-analysis-log_knee-ecog.pdf',width=10,height=3)







##### only DBS
db_fig_B02_06 <- db_fig_B02_04 %>% 
  filter(str_starts(electrode,'dbs_L'))

db_fig_B02_06_percent <- db_fig_B02_06 %>% 
  group_by(freq_range) %>% 
  summarise(f_offset = sum(abs(offset - offset_ref)<bw_offset)/length(offset),
            f_exponent = sum(abs(exponent - exponent_ref)<bw_exponent)/length(exponent),
            f_log_knee = sum(abs(log_knee - log_knee_ref)<bw_log_knee)/length(log_knee)) %>%
  mutate(label_offset=str_c(signif(f_offset*100,3),'%'),
         label_exponent=str_c(signif(f_exponent*100,3),'%'),
         label_log_knee=str_c(signif(f_log_knee*100,3),'%')) 

scale_fill_B02_06 <-
  scale_fill_gradientn(limits=c(0,100),oob=scales::squish,
                       colours=c('#DEEBF7','#C6DBEF','#9ECAE1','#6BAED6','#4292C6','#2171B5','#08519C','#08306B')) 

db_fig_B02_06 %>%
  ggplot() +
  aes(y=offset_ref,x=offset)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_offset,alpha=0.9) +
  scale_fill_B02_06 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-0.5,6),ylim=c(-0.5,6))+
  geom_text(aes(label=label_offset,x=NULL,y=NULL),x=-0.5,y=6,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_06_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_06a_senitivity-analysis-offset-dbs.pdf',width=10,height=3)


db_fig_B02_06 %>%
  ggplot() +
  aes(y=exponent_ref,x=exponent)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_exponent,alpha=0.9) +
  scale_fill_B02_06 +
  geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(0.5,5.5),ylim=c(0.5,5.5))+
  geom_text(aes(label=label_exponent,x=NULL,y=NULL),x=0.5,y=5.5,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_06_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_06b_senitivity-analysis-exponent-dbs.pdf',width=10,height=3)

db_fig_B02_06 %>%
  ggplot() +
  aes(y=log_knee_ref,x=log_knee)+
  geom_abline(slope=1,intercept = 0,color='gray60',size=2)+
  geom_hex(aes(fill = ..count..),binwidth=bw_log_knee,alpha=0.9) +
  scale_fill_B02_06 +
  #geom_density2d(breaks=c(0.05,0.25),adjust=c(1,1),color='#800000',alpha=0.5,size=0.5)+
  coord_equal(xlim=c(-1,2),ylim=c(-1,2))+
  geom_text(aes(label=label_log_knee,x=NULL,y=NULL),x=-1,y=2,vjust=1,hjust=0,color='gray60',
            data=db_fig_B02_06_percent,size=3)+
  facet_wrap(~freq_range,nrow=1)
ggsave('fig/B02_06c_senitivity-analysis-log_knee-dbs.pdf',width=10,height=3)



######

db_fig_B02_05 %>%
  filter(freq_range=='5-250hz') %>%
  mutate(has_knee = log_knee > 0, has_knee_ref = log_knee_ref > 0) %>%
  group_by(has_knee, has_knee_ref) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  mutate(f=n/sum(n))

#   has_knee has_knee_ref     n       f
# 1 FALSE    FALSE           83 0.00908
# 2 FALSE    TRUE            79 0.00864
# 3 TRUE     FALSE          591 0.0647 
# 4 TRUE     TRUE          8387 0.918 

#~7% modify their knee value more than 1.4dB 
#~7% gain a knee, 1% loses the knee. 


