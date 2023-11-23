library(tidyverse)
library(glue)

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

files <- dir(str_c(PATH_ANALISYS,'/data'),pattern=glob2rx('sub*-spectrum.csv'),full.names = TRUE)
names(files) <- basename(files)
spectra <- files %>% 
  map_df(read_csv, .id = "file_name")

spectra_l <- spectra %>%
  select(-'...1') %>%
  pivot_longer(cols=-c(file_name,freq),names_to='label',values_to='power') %>%
  mutate(subject=str_extract(file_name,'(EM\\d*)')) %>%
  mutate(run_str=str_extract(file_name,'run-(\\d*)')) %>%
  mutate(run_id=as.numeric(str_extract(run_str,'(\\d)+'))) %>%
  mutate(log_power = log10(power) + 12) %>% #volts^2 to uV^2
  left_join(elecloc) %>%
  filter(type %in% c('thal','ctx')) %>%
  filter(freq>0) %>%
  filter(freq<=500) %>%
  filter(!is.na(log_power)) 

spectra_l %>%
  ggplot() +
  aes(x=freq,y=log_power,color=type,group=label) +
  geom_line(alpha=0.2)+
  facet_wrap(~subject)+
  scale_x_log10()
ggsave('fig/A02_01_EP-spectra-500Hz.pdf')


LFP_spectra_agg_target <- spectra_l %>%
  mutate(sub_el = str_c(subject,'_',label)) %>%
  group_by(type,freq) %>%
  summarise(mean_power = mean(log_power, na.rm=TRUE), 
            sd_power = sd(log_power, na.rm=TRUE), 
            mad_power = mad(log_power, na.rm=TRUE), 
            median_power = median(log_power, na.rm=TRUE), 
            n=n(),
            n_sub = length(unique(subject)),
            n_el = length(unique(sub_el))) %>%
  mutate(sem=sd_power/sqrt(n)) %>%
  mutate(txt = str_c('EP[',type,'] Ns=',n_sub,' Ne=',n_el)) %>%
  write_tsv("fig/A02_02_LFP-group-level-sEEG-target.tsv")

color_ep_ctx = '#C4604F' 
color_ep_thal = '#36A5D1' #ET VIM

LFP_spectra_agg_target %>%
  ggplot() +
  aes(y=median_power,x=freq,ymax=median_power+mad_power,ymin=median_power-mad_power,color=type) +
  geom_line() +
  geom_smooth(stat="identity") +
  geom_text(aes(label=txt,x=NULL),x=0,y=-3,hjust=0,parse=FALSE,size=3,data=LFP_spectra_agg_target %>% filter(freq==2)) +
  scale_y_continuous(breaks=-4:3,labels=10^(-4:3)) +
  scale_x_log10() +
  facet_wrap(~type,ncol=1) +
  scale_color_manual(values=c(color_ep_ctx,color_ep_thal)) +
  theme(legend.position='none') +
  #coord_cartesian(ylim=c(-4,4),xlim=c(1,250))+
  labs(x='Frequency (Hz)', y="Power (uV^2/Hz)") +
  geom_vline(xintercept = c(250), linetype=2, color='gray')
ggsave("fig/A02_02_LFP-group-level-sEEG-target.pdf",width=4,height=5)






  


