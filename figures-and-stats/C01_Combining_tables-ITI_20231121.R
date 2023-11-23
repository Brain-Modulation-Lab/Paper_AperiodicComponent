library(tidyverse)
library(glue)
library(readr)
library(magrittr)

setwd("/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures")

PATH_DATA <- '/Volumes/Nexus/DBS'

# Produced phoneme
subjects<-
  tibble(subject=paste0('DBS',c(seq(3001,3032),seq(4057,4088)))) %>%
  mutate(path_annot=glue("{PATH_DATA}/{subject}/Preprocessed Data/Sync/annot/{subject}_inter_trial_interval.txt")) %>%
  mutate(exists_annot=file.exists(path_annot)) %>%
  filter(exists_annot)

annot_files <- subjects$path_annot
names(annot_files) <- subjects$subject
annot <- annot_files %>% 
  map_df(read_tsv, .id = "subject")


annot %>%
  filter(!is.na(duration)) %>%
  group_by(subject,session_id) %>%
  summarise(total_duration = sum(duration,na.rm=TRUE),n=n(),mean_duration=mean(duration)) %>%
  ungroup() %>%
  summarise(median_total_duration= median(total_duration), mad_tatal_duration = mad(total_duration),
            median_n= median(n), mad_n = mad(n),
            median_duration= median(mean_duration), mad_duration = mad(mean_duration))

# median_total_duration mad_tatal_duration median_n mad_n median_duration mad_duration
# <dbl>              <dbl>    <dbl> <dbl>           <dbl>        <dbl>
#   1                  254.               65.7      118  1.48            2.27        0.448


annot %>%
  group_by(subject,session_id) %>%
  summarise(total_duration = sum(duration,na.rm=TRUE)) %>%
  group_by(subject) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  summarise(median_n= median(n), mad_n = mad(n))
# median_n mad_n
# <int> <dbl>
#   1        3  1.48


