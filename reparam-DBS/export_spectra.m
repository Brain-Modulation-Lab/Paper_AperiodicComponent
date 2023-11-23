function export_spectra(SUBJECT)
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

%getting samplig rate
Fs = bml_getFs(D1);
%calculating welch window in seconds
% ww = min([produced_triplet.duration; stimulus_triplet.duration; ITI.duration]);
ww = MIN_SEGMENT_FOR_SPECTRUM;      

getspectrum = @(x,j) pwelch(x(j,:), Fs .* ww, [], [], Fs);

spectra_db=table();
fig = figure('Position',[0,0,1500,800]);

%% looping through sessions
for s=1:height(session)
  SESSION_ID=session.session_id(s);
  fprintf('subject %s session %i\n',SUBJECT,SESSION_ID);
  session_s = session(s,:);
  electrode_s = bml_annot_filter(electrode,session_s);

  %for each channel
  for l=1:length(D1.label)
    
    if ismember(D1.label{l},electrode_s.electrode)
      % Calculate a power spectrum with Welch's method

      cfg=[];
      cfg.channel = D1.label{l};
      cfg.trials = trial1.id(trial1.session_id==SESSION_ID);
      Dsl = ft_selectdata(cfg,D1);

      cfg=[];
      cfg.epoch = bml_annot_filter(trial,session_s);
      cfg.epoch = bml_annot_difference(cfg.epoch, artifact(strcmp(artifact.label, D1.label{l}),:));
      if height(cfg.epoch) < 10 %not enough trials after artifact rejection
        continue
      end
      cfg.epoch = cfg.epoch(cfg.epoch.duration >= MIN_SEGMENT_FOR_SPECTRUM,:);
      Dsl_cont = bml_redefinetrial(cfg,Dsl);      
      cfg=[];
      cfg.trials=find(~cellfun(@(x) any(isnan(x)), Dsl_cont.trial));
      Dsl_cont = ft_selectdata(cfg,Dsl_cont);
      if isempty(Dsl_cont.trial)
        continue
      end
      
      cfg=[];
      cfg.epoch = bml_annot_filter(produced_triplet,session_s);
      cfg.epoch = bml_annot_difference(cfg.epoch, artifact(strcmp(artifact.label, D1.label{l}),:));
      cfg.epoch = cfg.epoch(cfg.epoch.duration >= MIN_SEGMENT_FOR_SPECTRUM,:);
      Dsl_prod = bml_redefinetrial(cfg,Dsl);
      cfg=[];
      cfg.trials=find(~cellfun(@(x) any(isnan(x)), Dsl_prod.trial));
      Dsl_prod = ft_selectdata(cfg,Dsl_prod);
      
      cfg=[];
      cfg.epoch = bml_annot_filter(stimulus_triplet,session_s);
      cfg.epoch = bml_annot_difference(cfg.epoch, artifact(strcmp(artifact.label, D1.label{l}),:));
      cfg.epoch = cfg.epoch(cfg.epoch.duration >= MIN_SEGMENT_FOR_SPECTRUM,:);
      Dsl_stim = bml_redefinetrial(cfg,Dsl);
      cfg=[];
      cfg.trials=find(~cellfun(@(x) any(isnan(x)), Dsl_stim.trial));
      Dsl_stim = ft_selectdata(cfg,Dsl_stim);
      cfg=[];
      cfg.trials=find(cellfun(@(x) size(x,2) >= MIN_SEGMENT_FOR_SPECTRUM .* Fs, Dsl_stim.trial));
      Dsl_stim = ft_selectdata(cfg,Dsl_stim);
     
      cfg=[];
      cfg.epoch = bml_annot_filter(ITI,session_s);
      cfg.epoch = bml_annot_difference(cfg.epoch, artifact(strcmp(artifact.label, D1.label{l}),:));
      cfg.epoch = cfg.epoch(cfg.epoch.duration >= MIN_SEGMENT_FOR_SPECTRUM,:);
      Dsl_iti = bml_redefinetrial(cfg,Dsl);
      cfg=[];
      cfg.trials=find(~cellfun(@(x) any(isnan(x)), Dsl_iti.trial));
      Dsl_iti = ft_selectdata(cfg,Dsl_iti);

      %% fitting continuous data
      [~, freqs] = getspectrum(Dsl_cont.trial{1},1);
      
      Dslcw = bml_apply(getspectrum,Dsl_cont,1);
      Dslcw = [Dslcw.trial{:}];
      median_cont = median(Dslcw,2);
      tmp = cell2table(repmat({SUBJECT,SESSION_ID,D1.label{l},'cont'},length(freqs),1),...
           'VariableNames',{'subject','session_id','electrode','type'});
      tmp.frequency = freqs;
      tmp.power = median_cont;
      spectra_db = bml_annot_rowbind(spectra_db,tmp);

      %% fitting production epochs 
      Dslpw = bml_apply(getspectrum,Dsl_prod,1);
      Dslpw = [Dslpw.trial{:}];
      median_prod = median(Dslpw,2);
      %plot(freqs,log10(Dslpw),'Color',[hex2rgb('#888888') 0.1])
      tmp = cell2table(repmat({SUBJECT,SESSION_ID,D1.label{l},'prod'},length(freqs),1),...
           'VariableNames',{'subject','session_id','electrode','type'});
      tmp.frequency = freqs;
      tmp.power = median_prod;
      spectra_db = bml_annot_rowbind(spectra_db,tmp);
      
      %% fitting stimulus epochs
      Dslsw = bml_apply(getspectrum,Dsl_stim,1);
      Dslsw= [Dslsw.trial{:}];
      median_stim = median(Dslsw,2);
      %plot(freqs,log10(Dslsw),'Color',[hex2rgb('#888888') 0.1])
      tmp = cell2table(repmat({SUBJECT,SESSION_ID,D1.label{l},'stim'},length(freqs),1),...
           'VariableNames',{'subject','session_id','electrode','type'});
      tmp.frequency = freqs;
      tmp.power = median_stim;
      spectra_db = bml_annot_rowbind(spectra_db,tmp);
      
      
      %% fitting ITI epoch
      Dsliw = bml_apply(getspectrum,Dsl_iti,1);
      Dsliw= [Dsliw.trial{:}];
      median_iti = median(Dsliw,2);
      %plot(freqs,log10(Dsliw),'Color',[hex2rgb('#888888') 0.1])
      tmp = cell2table(repmat({SUBJECT,SESSION_ID,D1.label{l},'iti'},length(freqs),1),...
           'VariableNames',{'subject','session_id','electrode','type'});
      tmp.frequency = freqs;
      tmp.power = median_iti;
      spectra_db = bml_annot_rowbind(spectra_db,tmp);
      

      %% plotting
      clf(fig);

      % figure overlay of different epochs in log linear space
      subplot(1,2,1)
      plot(freqs,log10([median_prod,median_stim,median_iti]))
      xlim([2,250])
      ylim(YLIM)
      legend({'Prod','Stim','ITI'},'Location','northeast')
      hold on
      bands = bml_get_canonical_bands([2,250]);
      bands.fmid = 0.5.*(bands.fstarts + bands.fends);
      for i=1:height(bands)
        fill([bands.fstarts(i),bands.fstarts(i),bands.fends(i),bands.fends(i)],...
             YLIM(1) + [0,0.05.*range(YLIM),0.05.*range(YLIM),0],...
             hex2rgb(bands.color(i)),'EdgeColor','black','Marker','none');
        text(bands.fmid(i),YLIM(1) + 0.025.*range(YLIM),bands.symbol{i});
      end
      title([SUBJECT ' S' num2str(s) ' ' D.label{l}])

      % figure overlay of different epochs in log log space
      subplot(1,2,2)
      plot(log10(freqs),log10([median_prod,median_stim,median_iti]))
      xlim(log10([2,250]))
      ylim(YLIM)
      legend({'Prod','Stim','ITI'},'Location','northeast')
      hold on
      bands = bml_get_canonical_bands([2,250]);
      bands.fmid = sqrt(bands.fstarts .* bands.fends);
      for i=1:height(bands)
        fill(log10([bands.fstarts(i),bands.fstarts(i),bands.fends(i),bands.fends(i)]),...
             YLIM(1) + [0,0.05.*range(YLIM),0.05.*range(YLIM),0],...
             hex2rgb(bands.color(i)),'EdgeColor','black','Marker','none');
        text(log10(bands.fmid(i)),YLIM(1) + 0.025.*range(YLIM),bands.symbol{i});
      end
      xticks(log10(bands.fends))
      xticklabels(bands.fends)

      export_fig([PATH_FIG filesep SUBJECT '_S' num2str(s) '_' D.label{l}  '_fooof.png'],'-m3',fig);
      export_fig([PATH_FIG filesep SUBJECT '_S' num2str(s) '_' D.label{l}  '_fooof.png'],'-m3',fig);
    end

  end
end

writetable(spectra_db, [PATH_ANALYSIS '\data\' SUBJECT '_power_spectra.tsv'],'Delimiter','\t','FileType','Text');




