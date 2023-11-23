import os
import glob
import numpy as np
import pandas as pd
import mne

# Get EDF files to preprocess
edf_files = glob.glob('/Volumes/Nexus/sEEG/sub*/*.edf') + glob.glob('/Volumes/Nexus/sEEG/sub*/*.EDF')
# this file was cropped into two
edf_files.remove('/Volumes/Nexus/sEEG/sub-EM1002/EM1002_ses-EMU_task-resting_run-1.edf')

# only want to include channels from sEEG in further processing
df_localization = pd.read_csv('/Volumes/Nexus/Users/zouj/sEEG_data/electrode-localizations.tsv', delimiter='\t')

for file in edf_files:
    subj = os.path.basename(os.path.dirname(file)).replace('sub-', '')
    print('=== preprocessing ' + subj + " ===")
    
    # read raw file
    print('=== read raw data ===')
    # only want to process sEEG channels 
    sEEG_ch_list = list(df_localization.query("Subject == " + "'" + subj +"'")['label'])
    raw = mne.io.read_raw_edf(file, include = sEEG_ch_list, preload=True)

    # want to do common averaging on each individual electrode
    ch_lists = df_localization.query("Subject == " + "'" + subj +"'").\
        groupby('electrode')['label'].apply(list).tolist()

    CA_list = []
    for ch_list in ch_lists:
        print("=== CA re-referencing for one electrode ===")
        ch_list = [ch for ch in ch_list if ch in raw.info['ch_names']]
        CA, ref_data  = mne.set_eeg_reference(raw.copy().pick(ch_list), 'average', copy = True )
        CA_list.append(CA)

    print("=== combining raw files ===")
    CA_list[0].add_channels(CA_list[1:])

    # export file 
    output_file = os.path.join('/Volumes/Nexus/Users/zouj/sEEG_data/preprocessed', os.path.basename(file))
    print('=== exporting to ' + output_file + " ===")
    CA_list[0].export(output_file, overwrite = True)