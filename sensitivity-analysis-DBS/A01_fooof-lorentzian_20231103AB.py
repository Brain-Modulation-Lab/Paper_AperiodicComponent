import os
import glob
import numpy as np
import mne
import re
import pandas as pd
from fooof import FOOOFGroup

# Get EDF files of preprocessed data

psd = pd.read_table("/Volumes/Nexus/Users/busha/Analysis/2021-07-07-FOOOF-reparam/data/ITI_power_spectra.tsv",delimiter='\t')



for sub in psd.subject.unique():
    #sub = 'DBS3001'
    print('=== Fitting ' + sub + " with FOOOF ===")

    spectra = psd[psd['subject']==sub]

    print('=== Fit Spectra with FOOOF ===')
    fg = FOOOFGroup(peak_width_limits=[2,25],
                  max_n_peaks = 6,
                  min_peak_height = 0.15,
                  peak_threshold = 2,
                  aperiodic_mode = 'lorentzian',
                  regularization_weight = 0
                  )
    
    fg.fit(spectra.freqs, spectra.get_data(), freq_range = [1, 200])
    
    # exporting
    print('=== exporting FOOOFGroup ===')
    fg_output = os.path.join('/Volumes/Nexus/Users/zouj/sEEG_data/fooof_outputs/fooof_1-200hz_r0_bw1_pw25/fits',
                             os.path.basename(file).replace('-spectrum.h5', '-fooof_group.json'))
    
    fg.save(fg_output, save_results=True, save_settings=True, save_data=True)

    print('=== exporting plots ===')
    for ind in range(len(fg)):
        fooof = fg.get_fooof(ind)

        ch_name = spectra.info['ch_names'][ind]
        fooof_plot_path = os.path.join('/Volumes/Nexus/Users/zouj/sEEG_data/fooof_outputs/fooof_1-200hz_r0_bw1_pw25',
                                                   'plots',
                                                    os.path.basename(file).replace('-spectrum.h5', '-'+ch_name+'_report.png'))
        fooof.save_report(fooof_plot_path, plt_log = True)
        
