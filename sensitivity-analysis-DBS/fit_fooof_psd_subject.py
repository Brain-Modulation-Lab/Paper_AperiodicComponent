import os
import glob
import numpy as np
import mne
import re
import pandas as pd
import matplotlib.pyplot as plt
from fooof import FOOOFGroup

def fit_fooof_psd_subject(sub): #,PATH_ANALYSIS,freq_range = [1, 250]):    
    
    #PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-10-31-FOOOF-sensitivity-analysis/fooof_outputs/fooof_1-200hz_r0_bw1_pw25'
    #freq_range=[1, 200]
    
    #PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-10-31-FOOOF-sensitivity-analysis/fooof_outputs/fooof_1-150hz_r0_bw1_pw25'
    #freq_range=[1, 150]
    
    
    #PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-10-31-FOOOF-sensitivity-analysis/fooof_outputs/fooof_1-250hz_r0_bw1_pw25'
    #freq_range=[1, 250]
    
    #PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-10-31-FOOOF-sensitivity-analysis/fooof_outputs/fooof_30-250hz_r0_bw1_pw25'
    #freq_range=[30, 250]
    
    PATH_ANALYSIS = '/Users/ao622/Dropbox (Personal)/Lab-BML/Expan/2023-10-31-FOOOF-sensitivity-analysis/fooof_outputs/fooof_5-250hz_r0_bw1_pw25'
    freq_range=[5, 250]
    
    psd = pd.read_table("/Volumes/Nexus/Users/busha/Analysis/2021-07-07-FOOOF-reparam/data/ITI_power_spectra.tsv",delimiter='\t')
    freqs = psd.columns[6:psd.shape[1]].values.astype(float)    
    
    print('=== Fitting ' + sub + " with FOOOF ===")

    spectra = psd[psd['subject']==sub].reset_index()
    spectra_array = 10**spectra.iloc[:,6:psd.shape[1]].values

    print('=== Fit Spectra with FOOOF ===')
    fg = FOOOFGroup(peak_width_limits=[2,25],
                max_n_peaks = 6,
                min_peak_height = 0.15,
                peak_threshold = 2,
                aperiodic_mode = 'lorentzian',
                regularization_weight = 0
                )

    fg.fit(freqs, spectra_array, freq_range = freq_range)

    # exporting
    print('=== exporting FOOOFGroup ===')
    fg_output = os.path.join(PATH_ANALYSIS,'fits', sub + '-fooof_group.json')

    fg.save(fg_output, save_results=True, save_settings=True, save_data=True)

    df_aper = pd.DataFrame(np.array([np.append(np.append(fooof.aperiodic_params,fooof.r_squared),fooof.error) for fooof in fg.group_results]), columns=['offset', 'log_knee', 'exponent','R2','error'])
    df_aper = pd.concat([spectra.iloc[:,0:5],df_aper],axis=1,ignore_index=True)
    df_aper = df_aper.set_axis(['idx','subject','session_id','electrode','type','offset', 'log_knee', 'exponent','R2','error'],axis='columns')
    df_aper.to_csv(os.path.join(PATH_ANALYSIS,'fits', sub + '-fooof_aper.csv'))

    print('=== exporting plots ===')
    df_per = pd.DataFrame()
    for ind in range(len(fg)):
        fooof = fg.get_fooof(ind)

        if fooof.n_peaks_ > 0:
            df_per = df_per.append(
                    pd.concat([pd.concat([spectra.iloc[[ind],0:5]]*fooof.n_peaks_,axis=0,ignore_index=True),
                        pd.DataFrame(np.append(np.array(np.arange(1,fooof.peak_params_.shape[0]+1),ndmin=2).T,np.append(fooof.peak_params_,fooof.gaussian_params_,axis=1),axis=1), columns=['gauss_id','CF','PW','BW','mean','height','sd'])],
                    axis=1,ignore_index=True),
                ignore_index=True)
            
        session_id = spectra['session_id'][ind]
        ch_name = spectra['electrode'][ind]
        fooof_plot_path = os.path.join(PATH_ANALYSIS,'plots',sub+'_S'+ str(session_id) +'_'+ch_name+'_report.png')
        fooof.save_report(fooof_plot_path, plt_log = True)
        
    df_per = df_per.set_axis(['idx','subject','session_id','electrode','type','gauss_id','CF','PW','BW','mean','height','sd'],axis='columns')
    df_per.to_csv(os.path.join(PATH_ANALYSIS,'fits', sub + '-fooof_per.csv'))

