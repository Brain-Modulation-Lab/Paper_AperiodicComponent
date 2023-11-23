%% Loading paths
ft_defaults
bml_defaults
format long
clear

%% Configuration Variables and Paths
PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures';
PATH_DATA='/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2021-11-16-FOOOF-figures/data';
PATH_AVERAGE_MNI = '/Volumes/Nexus/DBS/DBS_subject_lists/MNI_ICBM_2009b_NLIN_ASYM/cortex/CortexLowRes_15000V.mat';
PATH_SUBCORT_ATLAS = '/Volumes/Nexus/Resources/STN-Atlas/atlas_index.mat';
PATH_SUBCORT_ATLAS_VIM = '/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/atlas_index.mat';


cd(PATH_ANALYSIS)
electrode = readtable('data/A01_DBS_aper_coord_dx.tsv','Delimiter', '\t', 'TreatAsEmpty', 'NA','FileType','text');

%loading cortical reconstructions
average_mni = load(PATH_AVERAGE_MNI);

subcort = load(PATH_SUBCORT_ATLAS);
% subcort_vim = load(PATH_SUBCORT_ATLAS_VIM);
% nii_vimi = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/91.nii.gz');
% nii_vime = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/94.nii.gz');
% nii_vimip = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/104.nii.gz');
% nii_vimep = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/DISTAL (Ewert 2017)/lh/122.nii.gz');
% 
% subcort_vimi_lh_fv = ea_nii2fv(nii_vimi);
% subcort_vime_lh_fv = ea_nii2fv(nii_vime);
% subcort_vimip_lh_fv = ea_nii2fv(nii_vimip);
% subcort_vimep_lh_fv = ea_nii2fv(nii_vimep);

%loading VL posterior ventral from Morel atlas
nii_vlpv = ea_load_nii('/Users/ao622/git/leaddbs/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/MorelAtlasICBM2009b (Jakab 2008)/lh/VLpv.nii.gz');
subcort_vlpv_lh_fv = ea_nii2fv(nii_vlpv);

color_et_ecog = '#C4604F';% #ET ECoG
color_pd_ecog = '#6F67A6';% #PD ECoG
color_ep_seeg = '#8A4F80';% #EP sEEG
color_pd_stn = '#F7924A';% #PD STN
color_pd_gpi = '#F9BD00';% #PD GPi
color_et_vim = '#36A5D1';% #ET VIM
color_ep_cm = '#9EB859';% #EP CM

%% plotting brain alone

figure;
hold on
patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on
view(-90,0)
axis off; axis equal
camlight('headlight','infinite');
print(gcf,[PATH_ANALYSIS '/fig/A02_01_mni-brain-lateral-view.png'],'-dpng','-r300')


%% ecog - PD

sel_electrode = electrode(strcmp(electrode.electrode_type,'ecog') & (electrode.mni_nonlinear_x < 0) & contains(electrode.dx,'PD') ,:);

figure;
hold on
patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x-15, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,'filled',...
  'MarkerFaceAlpha',1,'MarkerFaceColor',color_pd_ecog,'MarkerEdgeColor','k','LineWidth',0.01)
view(-90,0)
axis off; axis equal
camlight('headlight','infinite');
scalebar(0,70,-50, 10, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_02_mni_ecog_only_PD_Left.png'],'-dpng','-r300')


%% ecog - PD and ET


figure;
hold on
patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

sel_electrode = electrode(strcmp(electrode.electrode_type,'ecog') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dx,['PD']) ,:);
scatter3(sel_electrode.mni_nonlinear_x-14, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,'filled',...
  'MarkerFaceAlpha',1,'MarkerFaceColor',color_pd_ecog,'MarkerEdgeColor','k','LineWidth',0.1)

sel_electrode = electrode(strcmp(electrode.electrode_type,'ecog') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dx,['ET']) ,:);
scatter3(sel_electrode.mni_nonlinear_x-15, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,'filled',...
  'MarkerFaceAlpha',1,'MarkerFaceColor',color_et_ecog,'MarkerEdgeColor','k','LineWidth',0.01)


view(-90,0)
axis off; axis equal
camlight('headlight','infinite');
scalebar(0,70,-50, 10, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_03_mni_ecog_PD_and_ET_Left.png'],'-dpng','-r300')




%% dbs - STN

sel_electrode = electrode(strcmp(electrode.electrode_type,'dbs') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dbs_target,['STN']),:);

side=2; %left
figure;
hold on
patch('vertices', subcort.atlases.fv{1,side}.vertices, 'faces', subcort.atlases.fv{1,side}.faces,...
'FaceColor', [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,75,'filled','o',...
     'MarkerFaceColor',color_pd_stn,'MarkerEdgeColor','k','LineWidth',0.01)

view(-90,0)
axis off; axis equal
camlight('headlight');
scalebar(-10,-5,-13, 2, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_04_mni_DBS_STN.png'],'-dpng','-r300')

%% dbs - Gpi

sel_electrode = electrode(strcmp(electrode.electrode_type,'dbs') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dbs_target,['GPi']),:);

side=2; %left
figure;
hold on
patch('vertices', subcort.atlases.fv{5,side}.vertices, 'faces', subcort.atlases.fv{5,side}.faces,...
'FaceColor', [.6 .6 .7], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,75,'filled','o',...
     'MarkerFaceColor',color_pd_gpi,'MarkerEdgeColor','k','LineWidth',0.01)

view(-90,0)
axis off; axis equal
camlight('headlight');
scalebar(-10,-5,-13, 2, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_04b_mni_DBS_GPi.png'],'-dpng','-r300')



%% dbs - VIM

sel_electrode = electrode(strcmp(electrode.electrode_type,'dbs') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dbs_target,['VIM']),:);

side=2; %left
figure;
hold on

patch('vertices', subcort_vlpv_lh_fv.vertices, 'faces', subcort_vlpv_lh_fv.faces,...
'FaceColor', [.6 .7 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,75,'filled','o',...
     'MarkerFaceColor',color_et_vim ,'MarkerEdgeColor','k','LineWidth',0.01)

view(-90,0)
axis off; axis equal
camlight('headlight');
scalebar(-10,-5,-13, 2, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_05_mni_DBS_VIM.png'],'-dpng','-r300')





%% dbs - STN, GPi and VIM

%STN
sel_electrode = electrode(strcmp(electrode.electrode_type,'dbs') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dbs_target,['STN']),:);

side=2; %left
figure;
hold on
patch('vertices', subcort.atlases.fv{1,side}.vertices, 'faces', subcort.atlases.fv{1,side}.faces,...
'FaceColor',  [.7 .6 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,75,'filled','o',...
     'MarkerFaceColor',color_pd_stn,'MarkerEdgeColor','k','LineWidth',0.01)

   
%GPi

sel_electrode = electrode(strcmp(electrode.electrode_type,'dbs') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dbs_target,['GPi']),:);

side=2; %left
patch('vertices', subcort.atlases.fv{5,side}.vertices, 'faces', subcort.atlases.fv{5,side}.faces,...
'FaceColor', [.6 .6 .7], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,75,'filled','o',...
     'MarkerFaceColor',color_pd_gpi,'MarkerEdgeColor','k','LineWidth',0.01)

%VIM   
sel_electrode = electrode(strcmp(electrode.electrode_type,'dbs') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dbs_target,['VIM']),:);


patch('vertices', subcort_vlpv_lh_fv.vertices, 'faces', subcort_vlpv_lh_fv.faces,...
'FaceColor',  [.6 .7 .6], 'EdgeColor', 'none', 'FaceAlpha',0.5, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

scatter3(sel_electrode.mni_nonlinear_x, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,75,'filled','o',...
     'MarkerFaceColor',color_et_vim ,'MarkerEdgeColor','k','LineWidth',0.01)

   
view(-90,0)
axis off; axis equal
camlight('headlight');
scalebar(-10,-5,-13, 2, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_06_mni_DBS_STN_GPi_and_VIM.png'],'-dpng','-r300')


%==========

%% ecog - ET


figure;
hold on
patch('vertices', average_mni.Vertices, 'faces', average_mni.Faces,...
'FaceColor', [.9 .9 .9], 'EdgeColor', 'none', 'FaceAlpha',1, ...
'facelighting', 'gouraud', 'specularstrength', 0, 'ambientstrength', 0.5, 'diffusestrength', 0.5)
hold on

sel_electrode = electrode(strcmp(electrode.electrode_type,'ecog') & (electrode.mni_nonlinear_x < 0) & ismember(electrode.dx,['ET']) ,:);
scatter3(sel_electrode.mni_nonlinear_x-15, sel_electrode.mni_nonlinear_y, sel_electrode.mni_nonlinear_z,'filled',...
  'MarkerFaceAlpha',1,'MarkerFaceColor',color_et_ecog,'MarkerEdgeColor','k','LineWidth',0.01)

view(-90,0)
axis off; axis equal
camlight('headlight','infinite');
scalebar(0,70,-50, 10, 'mm')

print(gcf,[PATH_ANALYSIS '/fig/A02_07_mni_ecog__ET_Left.png'],'-dpng','-r300')

