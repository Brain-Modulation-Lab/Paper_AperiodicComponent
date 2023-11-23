SUBJECT='DBS3028';

%function export_spectra(SUBJECT)
%
% This protocol performs line noise filtering, rereferencing, defines trials and 
% exports spectra for fooof fittin

bml_defaults
set(0,'DefaultLegendAutoUpdate','off')

%
FOOOF_PEACK_WIDTH_LIMITS = [2, 100];
FOOOF_APERIODIC_MODE = 'knee';
FOOOF_MAX_N_PEAKS = 10;
FOOOF_F_RANGE = [2, 250];
MIN_SEGMENT_FOR_SPECTRUM = 1;  %s
CRITERIA = 'C';
YLIM = [-3,6];

%defining paths
%SUBJECT='DBS3001'; %used for development/debugging
PATH_DATA = 'Z:\DBS';
PATH_SUBJECT = [PATH_DATA filesep SUBJECT];
PATH_SYNC = [PATH_SUBJECT '\Preprocessed Data\Sync'];
PATH_PROTOCOL = 'Z:\Users\busha\Analysis\2021-07-07-FOOOF-reparam';

PATH_ANALYSIS = PATH_PROTOCOL;  
PATH_FIG = [PATH_ANALYSIS '\figures'];

DATE=datestr(now,'yyyymmdd');

fprintf('spectrum export Analysis for subject %s \n',SUBJECT)
  
cd(PATH_SYNC)
%loading annotation tables
session = bml_annot_read(['annot/' SUBJECT '_session.txt']);
%tf_epoch = bml_annot_read(['annot/' SUBJECT '_trial_epoch.txt']);
electrode = bml_annot_read(['annot/' SUBJECT '_electrode.txt']);
artifact = bml_annot_read(['annot/' SUBJECT '_artifact_criteria_' CRITERIA '.txt']);
empty_electrode = bml_annot_read(['annot/' SUBJECT '_empty_electrode.txt']);


%loading behavioral epochs
produced_triplet = bml_annot_read(['annot/' SUBJECT '_produced_triplet.txt']);

if isfile(['annot/' SUBJECT '_stimulus_triplet.txt'])
  stimulus_triplet = bml_annot_read(['annot/' SUBJECT '_stimulus_triplet.txt']);
else
  stimulus_triplet = produced_triplet(:,{'id','session_id','trial_id'});
  stimulus_triplet.starts = produced_triplet.starts - 2;
  stimulus_triplet.ends = produced_triplet.starts - 0.5;
end

if isfile(['annot/' SUBJECT '_inter_trial_interval.txt'])
  ITI = bml_annot_read(['annot/' SUBJECT '_inter_trial_interval.txt']);
  pre_ITI = bml_annot_rename(ITI(:,{'id','starts','ends','trial_id_post','session_id'}),'trial_id_post','trial_id');
  pre_ITI.starts_ITI = pre_ITI.starts;
  trial = bml_annot_left_join(produced_triplet(:,{'id','starts','ends','session_id','trial_id'}),pre_ITI);
  
	ITI = ITI(~ismissing(ITI.duration),:);
  ITI = ITI(ITI.duration>0,:);
else
	ITI = produced_triplet(:,{'id','session_id','trial_id'});
  ITI.starts = produced_triplet.starts - 4;
  ITI.ends = produced_triplet.starts - 2;
  stimulus_triplet.ends = produced_triplet.starts - 0.5;
  trial = produced_triplet;
  trial.starts_ITI = produced_triplet.starts - 4;
end
  
trial.starts = trial.starts_ITI;
trial = bml_annot_table(trial);
trial = trial(~ismissing(trial.duration),:);

%adjusting epochs
%extending production epoch 0.5 before onset
produced_triplet = bml_annot_extend(produced_triplet,0.5,0); 

cd(PATH_PROTOCOL)
  
%% loading continuous preprocessed data 
load([PATH_SUBJECT filesep 'Preprocessed Data' filesep 'FieldTrip' filesep SUBJECT '_ft_raw_session.mat'],'D');

%% rerefrencing

%combining artifacts with empty_electrode table
cfg=[];
cfg.groupby='label';
artifact_empty = bml_annot_union(cfg, artifact, empty_electrode);

cfg=[];
cfg.epoch=bml_annot_extend(trial,1,1);
[D1, trial1] = bml_redefinetrial(cfg,D);

%masking artifacts and empty_electrodes with NaNs
cfg=[];
cfg.annot=artifact_empty;
cfg.label_colname = 'label';
cfg.complete_trial = true; %masks entire trials
cfg.value=NaN;
D1 = bml_mask(cfg,D1);

%% Rereferencing

% macro bipolar reference to central
el_macro = electrode(electrode.type=="macro",:);
if ~isempty(el_macro)
  cfg=[];
  cfg.method = 'bipolar'; 
  cfg.label = unique(el_macro.electrode);
  cfg.refchan = {'macro_c'}; 
  cfg.refkeep = false; 
  D1 = bml_rereference(cfg,D1);
end

% dbs lead CAR 
el_dbs = electrode(electrode.type=="dbs",:);
if ~isempty(el_dbs)
  cfg=[];
  cfg.method = 'CAR'; 
  cfg.group = el_dbs.connector;
  cfg.label = el_dbs.electrode;
  D1 = bml_rereference(cfg,D1);
end

% ecog CTAR
el_ecog = electrode(electrode.type=="ecog",:);
if ~isempty(el_ecog)
  cfg=[];
  cfg.method = 'CTAR'; 
  cfg.percent = 50;
  cfg.group = el_ecog.connector;
  cfg.label = el_ecog.electrode;
  D1 = bml_rereference(cfg,D1);
end

%####

SESSION_ID = 4;
CH = {'ecog_248','ecog_249','ecog_250','dbs_L1','dbs_L2','dbs_L3'};
session_s = session(session.session_id==SESSION_ID,:);
cfg=[];
cfg.channel = CH;
cfg.trials = trial1.id(trial1.session_id==SESSION_ID);
Dsl = ft_selectdata(cfg,D1);

cfg=[];
cfg.epoch = bml_annot_filter(ITI,session_s);
cfg.epoch = cfg.epoch(cfg.epoch.duration >= MIN_SEGMENT_FOR_SPECTRUM,:);
Dsl_iti = bml_redefinetrial(cfg,Dsl);

cfg=[];
cfg.viewmode='vertical';
ft_databrowser(cfg,Dsl_iti)


