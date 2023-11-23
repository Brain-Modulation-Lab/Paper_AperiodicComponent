bml_defaults
ft_defaults

cd("/Volumes/Nexus/Users/zouj/sEEG_data/preprocessed")
cfg=[];
cfg.chantype=[];
cfg.pattern='*.edf';
r = bml_info_raw(cfg)

median(r.duration)
mad(r.duration,1)
%600 +/- 280 s
%10.0 +/- 4.6 min



%collapsing for each subject
