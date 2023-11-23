import os
import glob
import numpy as np
import mne

# Get EDF files of preprocessed data
edf_files = glob.glob('/Volumes/Nexus/Users/zouj/sEEG_data/preprocessed/*')

for file in edf_files:
    print('=== preprocessing ' + os.path.basename(file) + " ===")
    print('=== Read Data ===')
    signal = mne.io.read_raw_edf(file, preload=True)

    print('=== Calculate Spectra with Welchs ===')
    spectra = signal.compute_psd(method = 'welch',
                                 n_fft = 1024, # 1 second window
                                 n_overlap = 512 # 500 ms overlap
                                 )

    output_file = os.path.join('/Volumes/Nexus/Users/zouj/sEEG_data/spectra', os.path.basename(file)[:-4]+'-spectrum.h5' )
    print('=== exporting to ' + output_file + " ===")
    spectra.save(output_file)