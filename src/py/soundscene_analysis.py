# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 09:36:49 2023

@author: m_lot
"""

import os
import glob
import numpy as np
from scipy.signal import spectrogram
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from dsp.f_weight_T import A_weight_T
from utils import load_in
from dsp import noct

dark = False
if dark == True:
    plt.style.use('dark_background')

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['figure.autolayout'] = True

filepath = os.path.join(r"C:\Users\m_lot\OneDrive - University of Salford\Audio\Soundscapes\EigenScape\Lite-EigenScape")

dBAlevels = {"BusyStreet" : 80, "Park": 50, "QuietStreet": 65}


for file in glob.glob(filepath + "\\" + "*.flac"):
    
    cal_val = dBAlevels[os.path.basename(file.split(sep=".")[0])]

    x_Pa, fs = load_in.audio_caltnc(file, cal_val)
    
    xA_Pa = A_weight_T(x_Pa, fs)
    
    nperseg = 8192
    overlap = 0.5
    f, t, SxxA = spectrogram(xA_Pa[:, 0], fs, scaling='spectrum', nperseg=nperseg,
                             noverlap=overlap*nperseg)
    
    f_plt = noct.noctfnoct(63, 8000, 1)
    
    
    fig, ax = plt.subplots(figsize=(8, 6))
    plt.pcolormesh(t, f, 10*np.log10(SxxA/2e-5**2), cmap='viridis', shading='auto')
    plt.colorbar(label="dB(A) re 2e-5 Pa")
    plt.clim(cal_val-40, cal_val)
    plt.ylim([45,  11220])
    plt.yscale('log')
    plt.yticks(f_plt, labels=np.around(f_plt))
    plt.xticks(60*np.arange(0, 11), np.arange(0, 11))
    plt.xlim([0, 600 + 1/fs])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        left=False,      # ticks along the left edge are off
        right=False,         # ticks along the right edge are off
        labelleft=False) # labels along the left edge are off
    plt.xlabel("Elapsed time, mins")
    plt.ylabel("Frequency, Hz")
    plt.title(os.path.basename(file))
    
    if dark == True:
        plt.savefig(os.path.join(filepath,
                                 os.path.basename(file).replace(".", "").replace("flac", "")
                                 + "_n" + str(nperseg) + "_dk.png"),
                    dpi=600, format='png')
    else:
        plt.savefig(os.path.join(filepath,
                                 os.path.basename(file).replace(".", "").replace("flac", "")
                                 + "_n" + str(nperseg) + ".png"),
                    dpi=600, format='png')
    
    plt.show()
    