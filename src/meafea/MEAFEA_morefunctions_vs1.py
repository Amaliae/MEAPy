#!/usr/bin/env python
# coding: utf-8

# In[1]:


from statistics import mean
import json
import pickle
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import statannotations
from statannotations.Annotator import Annotator
from scipy.stats import wilcoxon
from scipy import stats
from scipy.stats import spearmanr, skew, boxcox, boxcox_normmax, norm, probplot, boxcox_normplot

from datetime import date

import seaborn as sns 
import os 
import re
from statsmodels.graphics.tsaplots import plot_acf

from scipy.stats import mannwhitneyu

from itertools import permutations, combinations
import scipy.signal
from scipy.stats import iqr
from statsmodels.stats.stattools import medcouple
from math import e 
get_ipython().run_line_magic('matplotlib', 'inline')
from itertools import combinations
from itertools import product
from itertools import chain
import scipy.stats as stats
import sys, importlib, os
import McsPy.McsData
import McsPy.McsCMOS
from McsPy import ureg, Q_
from lxml import etree
from matplotlib.patches import PathPatch
from sklearn.linear_model import LinearRegression
import statsmodels.api as slowess # to build a LOWESS model
from scipy.interpolate import interp1d 

from statsmodels.sandbox.stats.multicomp import multipletests
import statsmodels.api as sm
from statsmodels.tsa.api import ExponentialSmoothing,  SimpleExpSmoothing
import neo
import quantities as pq
from elephant import statistics

from dtw import dtw, accelerated_dtw
import matplotlib.patches as mpatches
from matplotlib import cm
import nbimporter

from pathlib import Path


from meafea.module_V5Radar_vs1 import *

from meafea.moduleV15_vs1 import *

from meafea.MEAFEA_functions_vs1 import *




from meafea.BurstDetection_lib_vs1 import *




import statsmodels.formula.api as smf

import random


from functools import reduce

from bokeh.io import output_file, save

import tsfel

from sklearn.feature_selection import VarianceThreshold
from sklearn.manifold import TSNE

from sklearn.decomposition import PCA

from sklearn.cluster import KMeans

from sklearn.metrics import silhouette_score
from yellowbrick.cluster import KElbowVisualizer

from sklearn.decomposition import PCA

import shutil


# In[2]:


def scatter_doses(plotparam, dfs, Param, X, hue,  hueorder, subhue, group_order,  figsize):
    '''df if experiment compound selected dataframe, by groupby'''
    
    
    #dfs=dfs.dropna()
    from statsmodels.formula.api import ols
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    
    compounds=np.unique(dfs[X].values)[0]
    experiment=np.unique(dfs['Experiment'].values)[0]
    
    hues=np.unique(dfs[hue])
    
    hueorder=[h for h in hueorder if h in hues]
    
    
    
    colordict=plotparam['colordict']
    
    if len(hueorder)>1:
        
        
        max_value=dfs[Param].mean()+5*(np.std(dfs[Param].values))

        plt.style.use('ggplot')

        Labels=hueorder

        groups=np.unique(dfs[subhue].values)

        combin=list(combinations(Labels, 2))
        
        print(combin, 'combinations')
        
        plt.rcParams.update({'font.size': plotparam['labelsize']})
        
        
        plt.rc('xtick',labelsize=18)
        plt.rc('ytick',labelsize=18)

        figures, axs=plt.subplots(len(groups), len(hue_order)-1, figsize=figsize,  sharey=plotparam['sharey'])
        
        axs = np.array([axs]).reshape(len(groups),len(hue_order)-1)
        
        dfall=pd.DataFrame()

        for cindex, comb in enumerate(zip(hue_order[1:], hue_order[:-1])):

            print(cindex, comb)

            if np.where(np.array(hueorder)==comb[0])[0]>np.where(np.array(hueorder)==comb[1])[0]:
                comb=comb[::-1]
                
                
            for indg, g in enumerate(group_order):
                
                dfnewgroup=pd.DataFrame()
                
                df=dfs[dfs[subhue]==g]
                
                max_value=df[Param].mean()+3*(np.std(df[Param].values))
                
                #min_value=df[Param].mean()-(np.std(df[Param].values))
                
                xdf=df[df[hue]==comb[0]]
                
               
                xch=xdf['Channel ID'].unique()
                
                ydf=df[df[hue]==comb[1]]
                
                ych=ydf['Channel ID'].unique()
                
                print(xdf.shape, ydf.shape)
                
                commonch=np.intersect1d(xch, ych)
                
                
                
                
                x=xdf[xdf['Channel ID'].isin(commonch)==True].sort_values(by='Channel ID')[Param].values
                
                y=ydf[ydf['Channel ID'].isin(commonch)==True].sort_values(by='Channel ID')[Param].values
                
                
                relativechange=(y-(x+0.001))/(x+0.001)
                
                absolutechange=y-x
                
                dfnewgroup['Relative Change' + ' '+Param]=relativechange
                
                dfnewgroup['Absolute Change'+' '+Param]=absolutechange
                
                dfnewgroup['Group']=g
                
                dfnewgroup['Channel ID']=commonch
                
                dfnewgroup['Dose Label']=str(comb[0])+' to '+str(comb[1])
                    
                dfall=pd.concat([dfall, dfnewgroup], axis=0)
             
                    
                    
                colorlist=np.repeat(colordict[g], len(x))
                axs[indg, cindex].scatter(x, y, color=colorlist, label=g)
                
                
                axs[indg, cindex].set_xlabel(comb[0], fontsize=plotparam['labelsize'])
                axs[indg, cindex].set_ylabel(comb[1],  fontsize=plotparam['labelsize'])
                #axs[indg, cindex].set_aspect('equal')
                
                model = LinearRegression(fit_intercept=False).fit(x.reshape(-1, 1), y)#sm.OLS(y, sm.add_constant(x)).fit()
                y_hat=model.predict(x.reshape(-1, 1))

                r2=r2_score(y, y_hat)  
                
                #axs[indg, cindex].set_title(Param+compounds+g, fontsize=20)
 
                axs[indg, cindex].plot(x, y_hat, color='grey')   ###label="r\u00b2 = {:.3f}".format(r2)
                axs[indg, cindex].plot(x, x, color='black', label='Unity', linestyle='dotted')
                axs[indg, cindex].set_xlim(-0.1,  max_value) 
                axs[indg, cindex].set_ylim(-0.1,  max_value)# Set x-axis limits
                axs[indg, cindex].legend(loc='upper left') 
                
                plt.tight_layout()
                
                
                
  
    
                
                      
        #plt.savefig(os.path.join(plotparam['savepath'], compounds+experiment+'Scatter'+'.png'), dpi=800,  bbox_inches='tight')
        
        plt.show()
        
        return dfall
                    


# In[3]:


def change_doses(df, Params, hue,  hueorder):
    '''df if experiment compound selected dataframe, by groupby'''
    
    

    
    hues=np.unique(df[hue])
    
    hueorder=[h for h in hueorder if h in hues]
    
   
    
    if len(hueorder)>1:
        
        Labels=hueorder


        combin=list(combinations(Labels, 2))
        
        print(combin, 'combinations')
        
        dfall=pd.DataFrame()

        for cindex, comb in enumerate(zip(hueorder[1:], hueorder[:-1])):

            print(cindex, comb)

            if np.where(np.array(hueorder)==comb[0])[0]>np.where(np.array(hueorder)==comb[1])[0]:
                
                comb=comb[::-1]
                dfnewgroup=pd.DataFrame()


                xdf=df[df[hue]==comb[0]]


                xch=xdf['Channel ID'].unique()

                ydf=df[df[hue]==comb[1]]

                ych=ydf['Channel ID'].unique()

                print(xdf.shape, ydf.shape)

                commonch=np.intersect1d(xch, ych)

                x=xdf[xdf['Channel ID'].isin(commonch)==True].sort_values(by='Channel ID')

                y=ydf[ydf['Channel ID'].isin(commonch)==True].sort_values(by='Channel ID')
                
                dfnewgroup['Channel ID']=x['Channel ID'].values

                dfnewgroup['Dose Label']=str(comb[0])+' to '+str(comb[1])

                for Param in Params:

                    xp=x[Param].values


                    yp=y[Param].values


                    relativechange=(yp-(xp+0.001))/(xp+0.001)

                    absolutechange=yp-xp

                    dfnewgroup['Relative Change' + ' '+Param]=relativechange

                    dfnewgroup['Absolute Change'+' '+Param]=absolutechange



                    dfall=pd.concat([dfall, dfnewgroup], axis=0)

    
        
        
        return dfall
                    


# In[4]:


def scatter_doses_well(plotparam, dfs, Param, X, hue,  hueorder, subhue, group_order,  figsize):
    '''df if experiment compound selected dataframe, by groupby'''
    
    
    #dfs=dfs.dropna()
    from statsmodels.formula.api import ols
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    
    compounds=np.unique(dfs[X].values)[0]
    experiment=np.unique(dfs['Experiment'].values)[0]
    
    hues=np.unique(dfs[hue])
    
    hueorder=[h for h in hueorder if h in hues]
    
    
    
    colordict=plotparam['colordict']
    
    if len(hueorder)>1:
        
        
        max_value=dfs[Param].mean()+5*(np.std(dfs[Param].values))

        plt.style.use('ggplot')

        Labels=hueorder

        groups=np.unique(dfs[subhue].values)

        combin=list(combinations(Labels, 2))
        
        print(combin, 'combinations')
        
        plt.rcParams.update({'font.size': plotparam['labelsize']})
        
        
        plt.rc('xtick',labelsize=18)
        plt.rc('ytick',labelsize=18)

        figures, axs=plt.subplots(len(groups), len(hue_order)-1, figsize=figsize,  sharey=plotparam['sharey'])
        
        axs = np.array([axs]).reshape(len(groups),len(hue_order)-1)
        
        dfall=pd.DataFrame()

        for cindex, comb in enumerate(zip(hue_order[1:], hue_order[:-1])):

            print(cindex, comb)

            if np.where(np.array(hueorder)==comb[0])[0]>np.where(np.array(hueorder)==comb[1])[0]:
                comb=comb[::-1]
                
                
            for indg, g in enumerate(group_order):
                
                dfnewgroup=pd.DataFrame()
                
                df=dfs[dfs[subhue]==g]
                
                max_value=df[Param].mean()+3*(np.std(df[Param].values))
                
                #min_value=df[Param].mean()-(np.std(df[Param].values))
                
                xdf=df[df[hue]==comb[0]]
                
               
                xch=xdf['Well ID'].unique()
                
                ydf=df[df[hue]==comb[1]]
                
                ych=ydf['Well ID'].unique()
                
                print(xdf.shape, ydf.shape)
                
                commonch=np.intersect1d(xch, ych)
                
                
                
                
                x=xdf[xdf['Well ID'].isin(commonch)==True].sort_values(by='Well ID')[Param].values
                
                y=ydf[ydf['Well ID'].isin(commonch)==True].sort_values(by='Well ID')[Param].values
                
                
                relativechange=(y-(x+0.001))/(x+0.001)
                
                absolutechange=y-x
                
                dfnewgroup['Relative Change' + ' '+Param]=relativechange
                
                dfnewgroup['Absolute Change'+' '+Param]=absolutechange
                
                dfnewgroup['Group']=g
                
                dfnewgroup['Well ID']=commonch
                
                dfnewgroup['Dose Label']=str(comb[0])+' to '+str(comb[1])
                    
                dfall=pd.concat([dfall, dfnewgroup], axis=0)
             
                    
                    
                colorlist=np.repeat(colordict[g], len(x))
                axs[indg, cindex].scatter(x, y, color=colorlist, label=g)
                
                
                axs[indg, cindex].set_xlabel(comb[0], fontsize=plotparam['labelsize'])
                axs[indg, cindex].set_ylabel(comb[1],  fontsize=plotparam['labelsize'])
                #axs[indg, cindex].set_aspect('equal')
                
                model = LinearRegression(fit_intercept=False).fit(x.reshape(-1, 1), y)#sm.OLS(y, sm.add_constant(x)).fit()
                y_hat=model.predict(x.reshape(-1, 1))

                r2=r2_score(y, y_hat)  
                
                #axs[indg, cindex].set_title(Param+compounds+g, fontsize=20)
 
                axs[indg, cindex].plot(x, y_hat, color='grey')   ###label="r\u00b2 = {:.3f}".format(r2)
                axs[indg, cindex].plot(x, x, color='black', label='Unity', linestyle='dotted')
                axs[indg, cindex].set_xlim(-0.1,  max_value) 
                axs[indg, cindex].set_ylim(-0.1,  max_value)# Set x-axis limits
                axs[indg, cindex].legend(loc='upper left') 
                
                plt.tight_layout()
                
                
                
  
    
                
                      
        #plt.savefig(os.path.join(plotparam['savepath'], compounds+experiment+'Scatter'+'.png'), dpi=800,  bbox_inches='tight')
        
        plt.show()
        
        return dfall
                    


# In[2]:


def get_file_and_chid(gdf, self):
    
    ###works with meaoobj
    
    
    ###cgdf is single row for a network burst
    
    experiment=gdf['Experiment'].unique()[0]
    dlabel=gdf['Dose Label'].unique()[0]
    welllabel=gdf['Well Label'].unique()[0]
    
    print(welllabel, 'welllabel')
    wellspikes=self.spikedata[(self.spikedata['Experiment']==experiment) & (self.spikedata['Dose Label']==dlabel)
                         & (self.spikedata['Well Label']==welllabel)]
    
    file=wellspikes['File'].unique()[0][:-10]
    full_times=self.stamp[file][str([dlabel])]
    
    dose_st=full_times['start']
    dose_end=full_times['stop']
    
    return file, wellspikes, int(dose_st), int(dose_end)


def time_reverse(x, start):
    
    x=(x-start)/50
    
    return int(x)
    
    
    

       


# In[3]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from statsmodels.robust import mad



def plot_signal_with_regions(signal, fs, high_regions, high_threshold):
    """
    Plot the signal with high-amplitude regions as shaded areas and thresholds.
    Parameters:
        signal (numpy.ndarray): Input signal.
        fs (float): Sampling frequency (Hz).
        high_regions (list): List of high-amplitude regions.
        high_threshold (float): High-amplitude threshold.
        low_threshold (float): Low-amplitude threshold.
    """
    time = np.arange(len(signal)) / fs

    plt.figure(figsize=(12, 6))
    plt.plot(time, signal, label="Signal", color="blue", linewidth=1)

    # Plot thresholds
    plt.axhline(y=high_threshold, color="green", linestyle="--", label="High Threshold")
    plt.axhline(y=-high_threshold, color="green", linestyle="--")
   
    # Highlight high-amplitude regions
    for i, region in enumerate(high_regions):
        start_time = region["Start Index"] / fs
        end_time = region["End Index"] / fs
        plt.axvspan(
            start_time,
            end_time,
            color=f"C{i % 10}",  # Cycle through 10 distinct colors
            alpha=0.3,
            label=f"Region {i + 1}" if i < 10 else None,
        )

    plt.title("Signal with High-Amplitude Regions")
    plt.xlabel("Time (seconds)")
    plt.ylabel("Amplitude")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.show()


# In[ ]:



# Function to normalize the signal
def normalize_signal(signal, original_fs, target_fs):
  
   
    # Normalize the signal
    normalized_signal = signal / np.max(np.abs(signal))
  
    num_samples = int(len(normalized_signal) * target_fs / original_fs)
    
    resampled_signal = resample(normalized_signal, num_samples)

    
    return resampled_signal 



# Function to find high amplitude regions
def find_high_amplitude_regions(signal, fs, threshold=0.5, chunk_size=100):
    
    
    """
    Identify regions in the signal where the amplitude exceeds the given threshold.
    
    
    """
    
    
    
    
    x_values = np.arange(len(signal))
    
    
    num_chunks = len(signal) // chunk_size + (1 if len(signal) % chunk_size else 0)
    
    
    high_amplitude_regions = []
    
    

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(signal))
        chunk_signal = signal[start_idx:end_idx]

        if np.max(np.abs(chunk_signal)) > threshold:
            marked_start = start_idx
            marked_end = end_idx
            marked_duration = (marked_end - marked_start) / fs
            marked_area = np.trapz(np.abs(chunk_signal), dx=1 / fs)
            max_abs_amplitude = np.max(np.abs(chunk_signal))
            high_amplitude_regions.append({
                "Start Index": marked_start,
                "End Index": marked_end,
                "Duration (seconds)": marked_duration,
                "Area Under Curve": marked_area,
                "Max Absolute Amplitude": max_abs_amplitude})

    return high_amplitude_regions


# In[6]:


def lowpass_filter(data, original_fs, cutoff_freq):
    
    from scipy.signal import butter, filtfilt
    """
    Apply a low-pass Butterworth filter to the data.

    Parameters:
        data (numpy.ndarray): The input signal.
        original_fs (float): Original sampling frequency (Hz).
        cutoff_freq (float): Low-pass filter cutoff frequency (Hz).

    Returns:
        numpy.ndarray: Low-pass filtered signal.
    """
    # Normalize the cutoff frequency to the Nyquist frequency
    nyquist_rate = original_fs / 2
    normalized_cutoff = cutoff_freq / nyquist_rate

    # Design a Butterworth low-pass filter
    b, a = butter(4, normalized_cutoff, btype='low', analog=False)

    # Apply the filter using filtfilt for zero-phase filtering
    filtered_data = filtfilt(b, a, data)
    return filtered_data

def downsample(data, downsampling_factor):
    """
    Downsample the data by selecting every nth sample.

    Parameters:
        data (numpy.ndarray): The input signal.
        downsampling_factor (int): Factor by which to downsample.

    Returns:
        numpy.ndarray: Downsampled signal.
    """
    # Downsample by selecting every nth sample
    return data[::downsampling_factor]

def process_lfp(data, original_fs, target_fs, lowpass_cutoff):
    """
    Low-pass filter and then downsample LFP data.

    Parameters:
        data (numpy.ndarray): The input LFP signal.
        original_fs (float): Original sampling frequency (Hz).
        target_fs (float): Desired downsampled frequency (Hz).
        lowpass_cutoff (float): Low-pass filter cutoff frequency (Hz).

    Returns:
        numpy.ndarray: Processed LFP signal.
    """
    # Calculate the downsampling factor
    downsampling_factor = int(original_fs / target_fs)
    if downsampling_factor < 1:
        raise ValueError("Target sampling rate must be lower than original sampling rate.")

    # Step 1: Apply low-pass filter
    filtered_data = lowpass_filter(data, original_fs, lowpass_cutoff)

    # Step 2: Downsample the filtered data
    downsampled_data = downsample(filtered_data, downsampling_factor)

    return downsampled_data


# In[4]:


def extract_raw(gdf, self, original_fs = 20000,
        target_fs = 2000,  
        lowpass_cutoff = 500, high_threshold=3,
                plotting=False, max_gap_high=0.15, 
                                                 max_gap_low=0.3, chunk_size=100):
    
    ###gdf is well and experiment and dose label groupped
    
    from statsmodels.robust import mad
    
    ##5 best bursts in one well
    
    
    sel_gdf=gdf.sort_values(by='N of Spike on Maximum', ascending=False)[:5]
    
    file, wellspikes, dose_st, dose_end=get_file_and_chid(gdf, self)
    
    Recordingfile=os.path.join(self.hdf5path, file)
        
    print(Recordingfile, 'Recordign File')
    
    recording=McsPy.McsData.RawData(Recordingfile)
    
    starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
    stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]
    esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
    
    dose_start=starts[np.where(esim==dose_st)[0]][0]
        
    dose_stop=stops[np.where(esim==dose_st)[0]][0]
    
    welldf=pd.DataFrame()
        
    for ix,  row in sel_gdf.iterrows():
        
       
        
    
        spikes=row['Spike on Maximum']
        
        spk_start=spikes[0]
        
        spk_end=spikes[-1]
        
        chlab=row['Maximum Electrode']
    
        chid=wellspikes[wellspikes['True Label']==chlab]['Channel ID'].unique()[0]
        
        index=int(recording.recordings[0].analog_streams[0].channel_infos[chid].row_index)
        
        v=recording.recordings[0].analog_streams[0].channel_data[index, dose_start:dose_stop]
        
        raw_start=time_reverse(spk_start, dose_st)
        
        raw_end=time_reverse(spk_end, dose_st)
        
        print(raw_start, raw_end)
        
        if raw_start>40000:
            
            substract=40000
            
        else:
            
            substract=raw_start
            
            
            
        
        burst=v[raw_start-substract:raw_start+40000].astype('float')
        
        if plotting==True:
            plt.figure()
            plt.plot(burst)
            plt.show()
        
        processed_data = process_lfp(burst, original_fs, target_fs, lowpass_cutoff)
        
        threshold=high_threshold*mad(processed_data)
        
        results=merge_regions_and_analyze_without_prominence(processed_data[20000:],
                                                             original_fs = original_fs,target_fs = target_fs,  
                                                             plotting=plotting,  high_threshold=high_threshold, 
                                                             max_gap_high=max_gap_high, 
                                                 max_gap_low=max_gap_low)
        
        res_df=pd.DataFrame([results])
        
        print(res_df.shape, 'resbefore')
        
        for col in ['Network_25th_Percentile_Spike_Time_Tiling_Coefficient',
                  'Network_Spike_Length_(ms)',
                  'Network_Mean_Spike_Time_Tiling_Coefficient',
                  'Duration on Maximum',
                  'MFR(Hz) on Maximum',
                  'Number_spikes_perelec_NetworkBurst_avg']:
            
             res_df[col]=row[col]
            
            
            
            
            
        
        
        
        
        print(res_df.shape, 'result_df')
        
        
        welldf=pd.concat([welldf, res_df], axis=0)
        
        print(welldf.shape, 'well_df')
        
        
        
        if plotting==True:
            
            plt.figure()

            plt.plot(processed_data)
            plt.axhline(y=threshold)
            plt.axhline(y=-1*threshold)
            
            plt.show()
        
        
    return welldf
       

dpath=r'G:\AH\iN47\Pharm'
meaobj=load_object(os.path.join(dpath,'ah301024detwithmad5fromcontrol_iN47_pharm.pkl'))
baselinedt=meaobj.spikedata[(meaobj.spikedata['Experiment']=='iN47CNVDelDIV48Rotenone1uMhr0') 
                           & (meaobj.spikedata['Dose Label']=='Control')]nets=meaobj.networkburstdata
nets['N of Spikes per channel']=nets['Spikes per channel'].apply(lambda x: [len(i) for i in x])
nets['Maximum Electrode']=nets[['Channel Labels', 'N of Spikes per channel']].apply(lambda x: x[0][np.argmax(x[1])], axis=1)
nets['Spike on Maximum']=nets[['Spikes per channel', 'N of Spikes per channel']].apply(lambda x: x[0][np.argmax(x[1])], axis=1)
nets['N of Spike on Maximum']=nets['Spike on Maximum'].apply(lambda x: len(x))netsmod=netsingle.groupby(['Experiment', 'Dose Label', 'Well Label', 'Group']).apply(lambda x: extract_raw(x, meaobj, 
                                                                                                     original_fs = 20000, 
        target_fs = 20000,  
        lowpass_cutoff = 10,  frequency_range1=[0, 50], 
            frequency_range=[0, 20], time_bandwidth=10, 
            num_tapers=None, window_params=[0.05, 0.01], 
                  min_nfft=0,
                  detrend_opt='off'))
# In[ ]:




cfg_file = tsfel.get_features_by_domain()
# In[7]:


def return_tdsfel_featuresinuse(usespectral, usestatistical, usetemporal, 
                                defaultSpectral=True, defaultStat=True, defaultTemp=True):
    
    cfg_file = tsfel.get_features_by_domain()
    
    cfg_filemod=cfg_file.copy()
    
    
    if defaultSpectral==True:

        usespectral=['FFT mean coefficient', 'Fundamental frequency',  'Max power spectrum', 'Maximum frequency',
                     'Median frequency',
                     'Power bandwidth', 'Spectral centroid', 
                     'Spectral entropy', 'Wavelet absolute mean', 
                     'Wavelet energy', 'Wavelet entropy']
        
    if defaultStat==True:
            usestatistical=list(cfg_file['statistical'].keys())
            
            
            
    if defaultTemp==True:
        
            usetemporal=['Area under the curve', 'Autocorrelation', 'Centroid', 'Mean absolute diff', 'Mean diff',
                     'Median absolute diff', 'Median diff', 'Negative turning points', 'Neighbourhood peaks', 
                     'Positive turning points', 'Signal distance', 'Slope', 'Sum absolute diff', 'Zero crossing rate']

    for feat in list(cfg_filemod['spectral'].keys()):

        if feat in usespectral:

            cfg_filemod['spectral'][feat]['use']='yes'

        else:

             cfg_filemod['spectral'][feat]['use']='no'

    for feat in list(cfg_filemod['statistical'].keys()):

        if feat in usestatistical:

            cfg_filemod['statistical'][feat]['use']='yes'

        else:

             cfg_filemod['statistical'][feat]['use']='no'

    for feat in list(cfg_filemod['temporal'].keys()):

        if feat in usetemporal:

            cfg_filemod['temporal'][feat]['use']='yes'

        else:

             cfg_filemod['temporal'][feat]['use']='no' 
        
    for feat in list(cfg_filemod['fractal'].keys()):


        cfg_filemod['fractal'][feat]['use']='no'

    for feat in list(cfg_filemod['fractal'].keys()):

        cfg_filemod['fractal'][feat]['use']='no'

    return cfg_filemod


# In[8]:


def Calculate_Spontant_Micro_rates(df, stamp, led, winsize=500, window=5000):
    
    dlabel=df['Dose Label'].unique()[0]

    file=df['File'].unique()[0]


    wellid=df['Well ID'].unique()[0]

    start, end=segment(df, stamp, dlabel, Duration=None)
    
    ledarg=check_led(file, led, stamp, df, 0, dlabel)
    
    arr=micro_rate(df, end, winsize)
    
    #print(ledarg, 'led arg')
    
    if ledarg==False:
        
        
        
        smoothed_df=smooth_micro_rate(arr, start, end, winsize=winsize, window=window)
        
        
    else:
        
        ledstart, ledend, leds, led_wells, light_label, light_dur=led_param(file, led, df, 
                                                                                        stamp, [], 0, dlabel)
        
        
        
        smoothed_df=smooth_micro_rate(arr, start, ledstart, winsize=winsize, window=window)
        
        
    return smoothed_df
        
        
        
def micro_rate(df, end, winsize):
    
    '''Calculates instanteous rates and the gradient'''
    
    
    ar=np.zeros(int((end-0)//winsize)+2)
    
    
    tmps=df['Timestamp [µs]'].values
    
    ####from spike 1 to spike 2 is rate of spike 1
    
    for tx, ty in zip(tmps[:-1], tmps[1:]):
        
        idxx=int(tx//winsize)
        
        idxy=int(ty//winsize)
        
        ##isi rate as ISI in 1 second
        
        rate=1000000/(ty-tx)
        
        
        ar[idxx:idxy]=rate
        
        
    return ar


def smooth_micro_rate(ar, start, end, winsize=500, window=5000):
    
    '''micro_rate can be further be smoothed with 2ms and 5ms windows'''
    


    ##time starts from zero   as   
    timeax=np.linspace(start, end+(2*winsize), len(ar[int(start//winsize):]))-start
    
    timeaxtrue=timeax+start
    

    dfroll=pd.DataFrame(ar[int(start//winsize):])
    
    timeroll=pd.DataFrame(timeax)
    ##Rolling MEAN (other stat when mean)
    #for each channel individually
    
    
    arr=dfroll.rolling(window//winsize, closed='both', axis=0).mean().values[int(window//winsize)-1:, :]
    
    timearr=timeroll.rolling(window//winsize, closed='both', axis=0).mean().values[int(window//winsize)-1:, :]
    #print(arr.shape, timearr.shape)
    
    ms5dfrolled=pd.DataFrame({'1/ISI': arr.reshape(arr.shape[0]), 'Time': timearr.reshape(timearr.shape[0])})
    
    ###m2df, 500us, 
    ms2df=pd.DataFrame({'1/ISI': ar[int(start//winsize):], 'Time': timeax, 'Time_from_start':timeaxtrue})
    
    
    return ms5dfrolled



def led_label(file, led, df, stamp, dlabel):
    
    
    
    pos=file.find('h5')+2
    
    
    
    key1=file[:pos]
    
    key2=dlabel
    
    dose_start, dose_end=segment(df, stamp, key2)
    
    if len(list(led[key1].keys()))>0:
        
        light_label=np.unique(np.array(led[key1]['Label']))
        
    else:
        
        light_label=''
    
    return light_label


def segment(df, stamp, Label, Duration=None):
    
    start=None
    end=None
    
    file=df['File'].unique()[0]
    
    pos=file.find('h5')+2
    
    
    key1=file[:pos]
    
    try: 
        start=int(stamp[key1][str([str(Label)])]['start'])
        end=int(stamp[key1][str([str(Label)])]['stop'])



        if Duration:

            end=start+(Duration*1000000)
            
    except:
        
        print('no file in stamp')
    
    return (start, end)
               
               
         


# In[10]:


def Calc_tsfel_features(ch, tsfelfeatures, fs, window_size=10):
    
    feat=pd.DataFrame()
    
    chdat=ch['1/ISI'].values
    
    if len(chdat)>10:
    
        feat=tsfel.time_series_features_extractor(tsfelfeatures, chdat, 
                                                  fs=fs, window_size=window_size)
    
    
    return feat

features_config=return_tdsfel_featuresinuse([],[], [], defaultSpectral=False)rate_data=meaobj.spikedata.groupby(meaobj.bychannel).apply(lambda x: Calculate_Spontant_Micro_rates(x, meaobj.stamp,
                                                                                        meaobj.LED, winsize=2000,
                                                                                        window=10000)).reset_index()

LED_dat['Time_Miliseconds']=LED_dat['Time_to_plot'].dt.total_seconds().astype(float)*1000

rate_data['Time_seconds']=rate_data['Time']/1000000

rtfeatures=rt.groupby(meaobj.bychannel).apply(lambda x: Calc_tsfel_features(x, features_config, 1, window_size=10)).reset_index()
# In[76]:




col=['Group']

col_order=['SBA', '22q11.2']

hue_order=['Control', '1µM 2h', '4µM 2h',
       '14µM 2h', '14µM 12h']
      

grouppalette={'SBA':'grey', '22q11.2':'dodgerblue'}



Set_Seaborn_Styles()
plotparam={}
plotparam['labelsize']=14
plotparam['savepath']=r'G:\\AH\\iN47\\Pharm\\Results'

plotparam['colordict']=grouppalette

plotparam['sharey']=False
# In[ ]:




#Params=['MFR(Hz)', 'BFR(Hz)', 'Median Burst_Duration_(sec)', 'Mean Burst_Duration_(sec)', 'Q1 Burst_Duration_(sec)', 
        'Max Burst_Duration_(sec)', 
        'Median Number_of_Spikes_per_Burst', 'Mean Number_of_Spikes_per_Burst', 'Q4 Burst_Duration_(sec)', 'Median First_to_Threshold_Spikes',
        'Q4 First_to_Threshold_Spikes', 'Q3 Burst_Min_Inter-Spike_Interval', 'Median Burstiness_Index', 'Median IBI']
# In[ ]:





# # Scatters control drug, pre post, this works for channels
# import matplotlib.patches as mpatches
legend_handles = [
    mpatches.Patch(color='grey', label='SBA'),
    mpatches.Patch(color='lightblue', label='22q11.2'),
    mpatches.Patch(edgecolor='black', facecolor='none', label='Rotenone', linestyle='-', linewidth=1),
    mpatches.Patch(edgecolor='black', facecolor='none', label='Vehicle', linestyle='--', linewidth=1),
]

# In[47]:


import matplotlib.ticker as ticker

for P in sorted_column_names_top_10:
    
    plt.figure()
    
    for d in CChall['Dose'].unique():
        
       # print(d)
        
    
        f=sns.lineplot(data=CChall[(CChall['Dose']==d) & (CChall['Dose Label'].isin(hue_order_relative)==True)
                                    & (CChall['Group'].isin(col_order)==True)
                                 & (CChall['Compounds']=='Vehicle')].dropna(),
                    y=P, x='Time',hue='Group', 
                hue_order=col_order, palette=grouppalette, markers='+', linestyle='--')

        f=sns.lineplot(data=CChall[(CChall['Dose']==d) & (CChall['Dose Label'].isin(hue_order_relative)==True) &
                                    (CChall['Group'].isin(col_order)==True)
                                 & (CChall['Compounds']=='Rotenone')].dropna(),
                    y=P, x='Time', hue='Group', 
                hue_order=col_order, palette=grouppalette, markers='o', linestyle='-')
   
    

    #f.set_xticklabels(hue_order_relative, fontsize=14,  ha='right', rotation_mode='anchor', rotation=45)
    
    plt.xticks(ticks=[ 0. ,  1. ,  2. ,  4. ,   5. ,  6. ,  8. ,    9. , 10 ,
       12 , 24], labels=[ 0 ,  1 ,  2 ,  4 ,    5 ,  6 ,  8 ,   9 , 10 ,
       12, 24])
    #f.xaxis.set_major_formatter(ticker.FixedFormatter(x[::5]))
    
    plt.legend(
    handles=legend_handles,
    title='Group',
    loc='center left',  # Position legend outside the plot
    bbox_to_anchor=(1, 0.5),  # Adjust position
    frameon=False,  # Add border around legend
    title_fontsize='10',
    fontsize='9', shadow=None
)for P in  sorted_column_names_top_10:
    
    aov = welch_anova(dv=P, between='Group',  data=dfChange2[dfChange2['Compounds']=='Rotenone'].dropna())
    print(aov, P)
    
    
    
for P in  sorted_column_names_top_10:
    
    print(P)
    
    print(pg.pairwise_gameshowell(data=dfChange2[dfChange2['Compounds']=='Rotenone'].dropna(), dv=P,
                        between='Group').round(3))
# In[25]:


def sort_describe_by_abs_mean_and_std(df, top_n=10):
    # Get the describe output
    describe_df = df.describe()
    
    # Sorting by max absolute mean value and min std
    sorted_columns = describe_df.loc[['50%', 'std']].T
    sorted_columns['abs_median'] = sorted_columns['50%'].abs()

    # Sort first by absolute mean (descending) and then by std (ascending)
    sorted_columns = sorted_columns.sort_values(by=['abs_median', 'std'], ascending=[False, True])

    # Return the top N column names
    return sorted_columns.head(top_n).index.tolist()


# Example: Apply this function to the previous DataFrame
#sorted_column_names_top_10 = sort_describe_by_abs_mean_and_std(dfChange3[dfChange3['Compounds']=='Rotenone'][use_params], top_n=10)
sorted_column_names_top_10=['Absolute Change Max Number_of_Spikes_per_Burst',
 'Absolute Change Q4 Number_of_Spikes_per_Burst', 'Absolute Change MFR(Hz)',
 'Relative Change MFR(Hz)', 'Relative Change BFR(Hz)', 
 'Relative Change Median Burstiness_Index', 
 'Absolute Change STD Number_of_Spikes_per_Burst',
 'Absolute Change Q3 Number_of_Spikes_per_Burst',
 'Absolute Change Mean Number_of_Spikes_per_Burst',
 'Relative Change Mean Number_of_Spikes_per_Burst', 
 'Relative Change Median Number_of_Spikes_per_Burst',
 'Relative Change Q1 First_to_Threshold_Spikes',
 'Relative Change Median First_to_Threshold_Spikes',
 'Relative Change Mean First_to_Threshold_Spikes', 
 'Relative Change Q3 Burst Duration(ms)' ]params_abs=['MFR(Hz)','STD Number_of_Spikes_per_Burst',
 'Q3 Number_of_Spikes_per_Burst',
 'Mean Number_of_Spikes_per_Burst',
 'Q1 First_to_Threshold_Spikes',
 'Median First_to_Threshold_Spikes',
'Q3 First_to_Threshold_Spikes', 'Median Burstiness_Index',
'BFR(Hz)',  'Q3 Burst_Duration_(sec)', 'Q3 Burst Duration(ms)']
# In[30]:


import pingouin as pg
from pingouin import welch_anova


# In[ ]:





# # Multivariate

# In[78]:


def Spider(burstpath, spikepath, joined, path, hue_order, plotby, compare, scaling):
    
    """This functon takes all params:
    1. MFR, BFR, Burst Paramms: for Control
    2. Change of MFR, BFR, Burst Params for Drugs
    And plots median of dist as radar plot
    
    Inputs arguments are: 
    Burst Folder path 
    Spike Folder path"""
    
    ##load
    
    if joined=='joined':
            
        dbfr=pd.read_csv(burstpath+'\\'+joined+'BFR(Hz).csv')
        dmfr=pd.read_csv(spikepath+'\\'+joined+'MFR(Hz).csv')
        dburst=pd.read_csv(burstpath+'\\'+joined+'Burst_analysis.csv')
        bfrchange=pd.read_csv(burstpath+'\\'+joined+'Change ofBFR(Hz).csv')
        mfrchange=pd.read_csv(spikepath+'\\'+joined+'Change ofMFR(Hz).csv') 

    else:
        
    
    
        dbfr=pd.read_excel(burstpath+'\\'+joined+'BFR(Hz).xlsx')
        dmfr=pd.read_excel(spikepath+'\\'+joined+'MFR(Hz).xlsx')
        dburst=pd.read_excel(burstpath+'\\'+joined+'Burst_analysis.xlsx')
        bfrchange=pd.read_csv(burstpath+'\\'+joined+'Change ofBFR(Hz).csv')
        mfrchange=pd.read_csv(spikepath+'\\'+joined+'Change ofMFR(Hz).csv') 

   
        
        
        
        
    
    
    

        
        
    
    ss=dburst.groupby(['Group', 'Experiment'])['Label'].value_counts().min()
    ndburst=dburst.groupby(plotby).apply(lambda x: x.iloc[random.sample(range(0, len(x)), ss), :]).reset_index(drop=True)
    for column in ['Duration', 'Spike Count', 'Max ISI', 'Min ISI', 'Mean ISI',
        'Burstiness', 'FT Spikes',
        'Surprise']:
    
    
        ndburst['Change'+column]=relative(ndburst, burstpath, hue_order, column, 'Group', 'Drugs', burstpath, 'False')['Change'+column]
    
    mergeon=['Experiment', 'Compound', 'Stim', 'Label',
       'Well ID', 'Group', 'Channel ID', 'LightWave', 'Age']
    rates=pd.merge(dbfr, dmfr, on=mergeon, how='inner')
    changerates=pd.merge(bfrchange, mfrchange, on=mergeon, how='inner')
    groupbys=['Experiment', 'Compound', 'Stim', 'LightWave', 'Group', 'Control']

    psetburst=['ChangeDuration', 'ChangeSpike Count',
       'ChangeMax ISI', 'ChangeMin ISI', 'ChangeMean ISI', 'ChangeBurstiness',
       'ChangeFT Spikes', 'ChangeSurprise']

    psetrates=['MFR(Hz)', 'BFR(Hz)']
    psetchrates=['ChangeMFR(Hz)', 'ChangeBFR(Hz)']
    
    ndburst.to_csv(path+'\\'+'ChangeofBurst_analysis.csv')
    
    
        
        
    
    if scaling=='True':
        
     
        scale=[-1, 1]



        for pset, df in zip([psetrates, psetchrates, psetburst, [p[6:] for p in psetburst]], [rates, changerates, ndburst, ndburst]):
    
            df['Control']=df['Label'].apply(lambda x: 'Yes' if x=='Control' else 'No')
    
            for p in pset:
        
                df[[p]]=df.groupby(groupbys)[[p]].transform(lambda x : np.apply_along_axis(MinMax, 0, x.values, scale))
    
    
    
    dflist=[]
    for pset, df in zip([psetchrates, psetburst], [changerates, ndburst]):
    
        df=df[pset+plotby].groupby(plotby).median().reset_index()
    
        dflist.append(df)
    
    dfChanges=pd.merge(dflist[0], dflist[1], on=plotby)
    
    
    dflist=[]

    psetburstc=[p[6:] for p in psetburst]

    for pset, df in zip([psetrates,  psetburstc], [rates,  ndburst]):
    
    
    
        df=df[df['Label']=='Control'][pset+plotby].groupby(plotby).median().reset_index()
    
        dflist.append(df)
    
    dfControl=pd.merge(dflist[0], dflist[1],  on=plotby)
    
    for dfp in [dfControl, dfChanges]:
    
    
        rrange=dfp.drop(plotby, axis=1).abs().max().values*1.5
        rangesxt=list(zip(rrange*-1, rrange))
        
        groupbie=[]
        
        for pb in plotby:
            
            if pb!=compare:
                
                groupbie.append(pb)
                
                
            
            
        
        dfp.groupby(groupbie).apply(lambda x: plot_radar(x, compare, groupbie, scaling, rangesxt, path))


# In[ ]:





# In[90]:


sorted_column_names_top_6=['Relative Change MFR(Hz)',
 'Relative Change BFR(Hz)',
 'Relative Change Median Burstiness_Index',
 
 'Relative Change Mean Number_of_Spikes_per_Burst',

 
 'Relative Change Q3 Burst Duration(ms)']

normalized_data = dt[dt['Compounds']=='Vehicle'].copy()
for param in sorted_column_names_top_6:
    normalized_data[param] = normalized_data[param].apply(lambda x: normalize(x, [normalized_data[param].min(), normalized_data[param].max()]))


# In[93]:


from scipy.stats import sem

radar_plot(normalized_data, sorted_column_names_top_6,  grouppalette, col_order, group_col='Group', rangeext=[-0.6, 0.6])
# In[86]:


def normalize(value, param_range):
    min_val, max_val = param_range
    return (value - min_val) / (max_val - min_val) * 2 - 1 
def round_angle(angle, base=60):
    return int(base * round(angle / base))


def radar_plot(data, Params, my_palette, group_order,  group_col='Group', rangeext=[-1.2, 1.2]):
    
    
    '''data used is relative change(*2)-1 to preserve direction '''
    
    # Radar plot
    categories = Params
    num_vars = len(categories)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]  # Complete the circle

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    
    group_stats = data.groupby(group_col)[Params].agg(['mean', sem])
    means = group_stats.xs('mean', level=1, axis=1)
    errors = group_stats.xs('sem', level=1, axis=1).multiply(1.96)  # 95% CI


     # Plot each group
    for group in group_order:
        
        values = means.loc[group, categories].tolist()
        values += values[:1]  # Complete the circle
        ci_lower = (means.loc[group, categories] - errors.loc[group, categories]).tolist()
        ci_upper = (means.loc[group, categories] + errors.loc[group, categories]).tolist()
        ci_lower += ci_lower[:1]
        ci_upper += ci_upper[:1]

        # Plot mean line
        ax.plot(angles, values, label=group, color=my_palette[group])
        # Fill area for confidence interval
        ax.fill_between(angles, ci_lower, ci_upper, alpha=0.2, color=my_palette[group])

    ax.set_yticks([rangeext[0], rangeext[0]/2, 0, rangeext[1], rangeext[1]/2])
    ax.set_xticks(angles[:-1])
    
    

    annotation_text = "\n".join(
        [f"{round_angle(np.degrees(angle))}° - {category}" for angle, category in zip(angles[:-1], categories)]
    )

    # Add text box to the plot
    plt.text(
        1.1, 0.5, annotation_text,
        transform=ax.transAxes, fontsize=14,
        verticalalignment='center', horizontalalignment='left',
        bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.5')
    )
    ax.set_ylim(rangeext)
    #ax.set_ylim(-1, 1)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

   


# In[140]:


def normalize(value, param_range):
    min_val, max_val = param_range
    return (value - min_val) / (max_val - min_val) * 2 - 1 


# In[144]:


def plot_radar(dt, compare, Params, Norm, rangesext, my_palette):
    
    gp=np.unique(dt[compare].values)
    
    columns=dt.columns
    
    #print(columns)
    
    variables=Params
    
    if Norm=='False':
        
        suf='Norm'
        
        ranges=[rangesext for v in variables]
    else:
        
        suf='Sym'
        
        ranges=[]
        
        for variable in variables:
            
            rmin= dt[variable].min()
            rmax=dt[variable].max()
            
            rangev=[dt[variable].min(),dt[variable].max()]
            
            ranges.append(rangev)
       
        
        
    #print(len(variables), len(ranges))
    
    #print(dlabel)
    
    fig1 = plt.figure(figsize=(6, 6))
    
    #plt.title(dlabel)
                      
                      
    
    #my_palette = plt.cm.get_cmap("PiYG", len(dt.index))
    
    handles=[]
    
    for ind, g in enumerate(gp): 
        
        data = dt[dt[compare]==g][Params]
        
        normalized_data = data.copy()
        
        for param in variables:
            
            normalized_data[param] = data[param].apply(lambda x: normalize(x, [data[param].min(),data[param].max()]))

        
        ##print(data, variables.values, ranges, my_palette(ind))
       
        radar = ComplexRadar(fig1, variables, ranges)
        radar.plot(normalized_data[variables].median().values, g,  color=my_palette[g], alpha=1)
        #radar.fill(data, color=my_palette(ind), alpha=0.4)
        r_patch = mpatches.Patch(color=my_palette[g], label=g)
        handles.append(r_patch)
        
    plt.legend(handles=handles, loc='best', bbox_to_anchor=(0.1, 0))
    
   

    #fig1.savefig(tosub+'/'+d+'Radarby'+compare+suf+'.svg', format='svg', bbox_inches='tight')
    plt.show()    
        


# In[9]:


def custom_legend_withN(df, color_map, cutomn=0):
    
    import matplotlib.patches as mpatches
    
    
    patches = []
    unique_groups=df["Group"].unique()
    group_counts = df["Group"].value_counts()  # sample sizes by group
    for group in unique_groups:
        group_color = color_map[group]
        if cutomn>0:
            n=cutomn
            
        else:
            n = group_counts[group]
        label_str = f"{group} (n={n})"
        patch = mpatches.Patch(color=group_color, label=label_str)
        patches.append(patch)

    return patches


# In[134]:


class ComplexRadar():
    def __init__(self, fig, variables, ranges,
                 n_ordinate_levels=4):
        
        angles = np.arange(0, 360, 360./len(variables))

        axes = [fig.add_axes([0.1,0.1,0.9,0.9], polar=True,
                label = "axes{}".format(i))
    
                for i in range(len(variables))]
        
        
        l, text = axes[0].set_thetagrids(angles, 
                                         labels=[v[15:] if 'Change' in v else v for v in variables], fontsize=12)
        [txt.set_rotation(angle-90) for txt, angle 
             in zip(text, angles)]
        axes[0].tick_params(axis='both', which='major', pad=20)
        for ax in axes[1:]:
            ax.patch.set_visible(False)
            ax.grid("off")
            ax.xaxis.set_visible(False)
        for i, ax in enumerate(axes):
            grid = np.linspace(*ranges[i], 
                               num=n_ordinate_levels)
            gridlabel = ["{}".format(round(x,2)) 
                         for x in grid]
            if ranges[i][0] > ranges[i][1]:
                grid = grid[::-1] # hack to invert grid
                          # gridlabels aren't reversed
            gridlabel[0] = "" # clean up origin
            ax.set_rgrids(grid, labels=gridlabel,
                         angle=angles[i])
            #ax.spines["polar"].set_visible(False)
            ax.set_ylim(*ranges[i])
        # variables for plotting
        
        ##print(angles, angles.shape, np.deg2rad(np.r_[angles, angles[0]], )
        self.angle = np.deg2rad(np.r_[angles, angles[0]])
        self.ranges = ranges
        self.ax = axes[0]
        
    def plot(self, data, g,  *args, **kw):
        
       
        ##print(self.angle, data, sdata[0])
        self.ax.plot(self.angle, np.r_[data, data[0]], label=g,  *args, **kw)
    def fill(self, data, *args, **kw):
        #sdata = _scale_data(data, self.ranges)
        self.ax.fill(self.angle, np.r_[data, data[0]], *args, **kw)


# In[99]:


def _scale_data(data, ranges):
    """scales data[1:] to ranges[0],
    inverts if the scale is reversed"""
    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        assert (y1 <= d <= y2) or (y2 <= d <= y1)
    x1, x2 = ranges[0]
    d = data[0]
    if x1 > x2:
        d = _invert(d, (x1, x2))
        x1, x2 = x2, x1
    sdata = [d]
    for d, (y1, y2) in zip(data[1:], ranges[1:]):
        if y1 > y2:
            d = _invert(d, (y1, y2))
            y1, y2 = y2, y1
        sdata.append((d-y1) / (y2-y1) 
                     * (x2 - x1) + x1)
    return sdata


# In[11]:


d


# In[12]:






def plot_signal_with_regions(signal, fs, high_regions, high_threshold):
    """
    Plot the signal with high-amplitude regions as shaded areas and thresholds.
    Parameters:
        signal (numpy.ndarray): Input signal.
        fs (float): Sampling frequency (Hz).
        high_regions (list): List of high-amplitude regions.
        high_threshold (float): High-amplitude threshold.
        low_threshold (float): Low-amplitude threshold.
    """
    time = np.arange(len(signal)) / fs

    plt.figure(figsize=(12, 6))
    plt.plot(time, signal, label="Signal", color="blue", linewidth=1)

    # Plot thresholds
    plt.axhline(y=high_threshold, color="green", linestyle="--", label="High Threshold")
    plt.axhline(y=-high_threshold, color="green", linestyle="--")
   
    # Highlight high-amplitude regions
    for i, region in enumerate(high_regions):
        start_time = region["Start Index"] / fs
        end_time = region["End Index"] / fs
        plt.axvspan(
            start_time,
            end_time,
            color=f"C{i % 10}",  # Cycle through 10 distinct colors
            alpha=0.3,
            label=f"Region {i + 1}" if i < 10 else None,
        )

    plt.title("Signal with High-Amplitude Regions")
    plt.xlabel("Time (seconds)")
    plt.ylabel("Amplitude")
    plt.legend(loc="upper right")
    plt.grid(True)
    plt.show()


# In[13]:





# In[15]:


def lowpass_filter(data, original_fs, cutoff_freq):
    from scipy.signal import butter, filtfilt
    """
    Apply a low-pass Butterworth filter to the data.

    Parameters:
        data (numpy.ndarray): The input signal.
        original_fs (float): Original sampling frequency (Hz).
        cutoff_freq (float): Low-pass filter cutoff frequency (Hz).

    Returns:
        numpy.ndarray: Low-pass filtered signal.
    """
    # Normalize the cutoff frequency to the Nyquist frequency
    nyquist_rate = original_fs / 2
    normalized_cutoff = cutoff_freq / nyquist_rate

    # Design a Butterworth low-pass filter
    b, a = butter(4, normalized_cutoff, btype='low', analog=False)

    # Apply the filter using filtfilt for zero-phase filtering
    filtered_data = filtfilt(b, a, data)
    return filtered_data

def downsample(data, downsampling_factor):
    """
    Downsample the data by selecting every nth sample.

    Parameters:
        data (numpy.ndarray): The input signal.
        downsampling_factor (int): Factor by which to downsample.

    Returns:
        numpy.ndarray: Downsampled signal.
    """
    # Downsample by selecting every nth sample
    return data[::downsampling_factor]

def process_lfp(data, original_fs, target_fs, lowpass_cutoff):
    """
    Low-pass filter and then downsample LFP data.

    Parameters:
        data (numpy.ndarray): The input LFP signal.
        original_fs (float): Original sampling frequency (Hz).
        target_fs (float): Desired downsampled frequency (Hz).
        lowpass_cutoff (float): Low-pass filter cutoff frequency (Hz).

    Returns:
        numpy.ndarray: Processed LFP signal.
    """
    # Calculate the downsampling factor
    downsampling_factor = int(original_fs / target_fs)
    if downsampling_factor < 1:
        raise ValueError("Target sampling rate must be lower than original sampling rate.")

    # Step 1: Apply low-pass filter
    filtered_data = lowpass_filter(data, original_fs, lowpass_cutoff)

    # Step 2: Downsample the filtered data
    downsampled_data = downsample(filtered_data, downsampling_factor)

    return downsampled_data


# In[16]:



# Function to normalize the signal
def normalize_signal(signal, original_fs, target_fs):
  
   
    # Normalize the signal
    normalized_signal = signal / np.max(np.abs(signal))
  
    num_samples = int(len(normalized_signal) * target_fs / original_fs)
    
    resampled_signal = resample(normalized_signal, num_samples)

    
    return resampled_signal 



# Function to find high amplitude regions
def find_high_amplitude_regions(signal, fs, threshold=0.5, chunk_size=100):
    
    
    """
    Identify regions in the signal where the amplitude exceeds the given threshold.
    
    
    """
    
    
    
    
    x_values = np.arange(len(signal))
    
    
    num_chunks = len(signal) // chunk_size + (1 if len(signal) % chunk_size else 0)
    
    
    high_amplitude_regions = []
    
    

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(signal))
        chunk_signal = signal[start_idx:end_idx]

        if np.max(np.abs(chunk_signal)) > threshold:
            marked_start = start_idx
            marked_end = end_idx
            marked_duration = (marked_end - marked_start) / fs
            marked_area = np.trapz(np.abs(chunk_signal), dx=1 / fs)
            max_abs_amplitude = np.max(np.abs(chunk_signal))
            high_amplitude_regions.append({
                "Start Index": marked_start,
                "End Index": marked_end,
                "Duration (seconds)": marked_duration,
                "Area Under Curve": marked_area,
                "Max Absolute Amplitude": max_abs_amplitude})

    return high_amplitude_regions


# In[3]:


def merge_regions_and_analyze_without_prominence(signal, original_fs=20000, target_fs=20000, 
                                                 high_threshold=0.5, low_threshold=0.15, 
                                                 max_gap_high=0.15, 
                                                 max_gap_low=0.3, chunk_size=100, 
                                                 plotting=True):
    """
    Process signal to identify high-amplitude regions, merge regions as per thresholds and gaps, and analyze the result,
    focusing on peaks with the highest positive and negative amplitudes without using prominence.
    """
    # Identify high-amplitude regions
    
    
    from scipy.signal import find_peaks
    
    from statsmodels.robust import mad


    from scipy.signal import resample
    
    fs=target_fs
    
    #norm_signal=normalize_signal(signal, original_fs, target_fs)
    
    threshold=mad(signal[:10000])*high_threshold
    
    low_threshold=mad(signal[:10000])*high_threshold/1.5
    
    
    
    
    
    high_amplitude_regions = find_high_amplitude_regions(signal, target_fs, 
                                                         threshold=threshold, chunk_size=chunk_size)
    
   

    
    plot_signal_with_regions(signal, fs, high_amplitude_regions, threshold)

    # Merge high-amplitude regions within the max gap
    merged_high_regions = []
    current_region = None

    for region in high_amplitude_regions:
        if current_region is None:
            current_region = region
        elif region["Start Index"] <= current_region["End Index"] + int(max_gap_high * fs):
            current_region["End Index"] = max(current_region["End Index"], region["End Index"])
            current_region["Duration (seconds)"] = (
                (current_region["End Index"] - current_region["Start Index"]) / fs
            )
            current_region["Area Under Curve"] += region["Area Under Curve"]
            current_region["Max Absolute Amplitude"] = max(
                current_region["Max Absolute Amplitude"], region["Max Absolute Amplitude"]
            )
        else:
            merged_high_regions.append(current_region)
            current_region = region
    if current_region:
        merged_high_regions.append(current_region)

    # Select the region with the maximum duration
    if not merged_high_regions:
        return None
    max_duration_high_region = max(merged_high_regions, key=lambda r: r["Duration (seconds)"])

    # Merge subsequent low-amplitude regions within the max gap
    combined_region = max_duration_high_region.copy()
    
    low_amplitude_regions=find_high_amplitude_regions(signal, fs, threshold=low_threshold, chunk_size=chunk_size)
    
    plot_signal_with_regions(signal, fs, low_amplitude_regions, threshold)
    
    for region in low_amplitude_regions:
        if (
            region["Start Index"] > combined_region["End Index"]
            and region["Start Index"] <= combined_region["End Index"] + int(max_gap_low * fs)
        ):
            combined_region["End Index"] = max(combined_region["End Index"], region["End Index"])
            combined_region["Duration (seconds)"] = (
                (combined_region["End Index"] - combined_region["Start Index"]) / fs
            )
            combined_region["Area Under Curve"] += region["Area Under Curve"]

    # Analyze parameters
    region_signal = signal[combined_region["Start Index"]:combined_region["End Index"]]
    
    
    ###find all peaks
    # Find all positive peaks above the threshold
    positive_peaks_indices, _ = find_peaks(region_signal, height=low_threshold)
    positive_peaks = region_signal[positive_peaks_indices]

    # Find all negative peaks below the negative threshold by inverting the signal
    negative_peaks_indices, properties = find_peaks(-region_signal, height=low_threshold)
    negative_peaks = region_signal[negative_peaks_indices]
    
    # Filter negative peaks to ensure they are below -threshold
    negative_peaks = negative_peaks[negative_peaks < -low_threshold]
    negative_peaks_indices = negative_peaks_indices[:len(negative_peaks)] 
    
    # Adjust indices accordingly
    
    negindeces=negative_peaks.argsort()
    
    ##sort from high to low
    
    negative_peaks_indicesperorder=negative_peaks_indices[negindeces]
    
    
    
    

    # Identify first positive and first negative peaks
    first_positive_peak_idx = positive_peaks_indices[0] if len(positive_peaks_indices) > 0 else None
    first_negative_peak_idx = negative_peaks_indices[0] if len(negative_peaks_indices) > 0 else None

    # Identify last negative peak
    last_negative_peak_idx = negative_peaks_indices[-1] if len(negative_peaks_indices) > 0 else None

    # Find maximum positive and negative peaks
    max_positive_peak_idx = np.argmax(region_signal)
    max_negative_peak_idx = np.argmin(region_signal)
    max_positive_peak = region_signal[max_positive_peak_idx]
    max_negative_peak = region_signal[max_negative_peak_idx]
    
    ##find a baseline.
    
    baseline_idx=max_positive_peak_idx+combined_region["Start Index"]
    
    trace=signal[int(baseline_idx-(1*fs)):int(baseline_idx+(1*fs))]
    
    baselinemedian=np.median(signal[np.where(signal[:baseline_idx]<low_threshold)[0]])
    
    
    
    ##afterpeaknegativepeakindices
    second_negative_peak_idx=negative_peaks_indicesperorder[negative_peaks_indicesperorder> max_positive_peak_idx][0] if len(negative_peaks_indicesperorder[negative_peaks_indicesperorder> max_positive_peak_idx])>0 else None
    
    first_negative_peak_idx=negative_peaks_indicesperorder[negative_peaks_indicesperorder< max_positive_peak_idx][0] if len(negative_peaks_indicesperorder[negative_peaks_indicesperorder< max_positive_peak_idx])>0 else None
    
    if second_negative_peak_idx!=None:
    
        aftersecindexes=np.where(abs(region_signal[second_negative_peak_idx:])
                                         <low_threshold)[0]+second_negative_peak_idx
        
        if len(aftersecindexes)>0:

            aftersec_idx=aftersecindexes.min()

            aftersecbaseline=np.median(region_signal[aftersecindexes.tolist()])
            
        else:
            second_negative_peak_idx=None
    
    if max_negative_peak_idx!=None:
        
        aftermaxindexes=np.where(abs(region_signal[max_negative_peak_idx:])
                                         <low_threshold)[0]+max_negative_peak_idx
        
        if len(aftermaxindexes)>0:

            aftermax_idx=aftermaxindexes.min()
            aftermaxbaseline=np.median(region_signal[aftermaxindexes.tolist()])
        
        else:
            max_negative_peak_idx=None
            aftermaxbaseline=None

       
    
    

    # Calculate rising rate
    if first_negative_peak_idx is not None and max_positive_peak_idx is not None and first_negative_peak_idx < max_positive_peak_idx:
        # Rising rate between first negative and first positive peaks
        ffrising_rate = (
            region_signal[max_positive_peak_idx] - region_signal[first_negative_peak_idx]
        ) / ((max_positive_peak_idx - first_negative_peak_idx) / fs)
        
        ffmaxtofirstratio=abs(region_signal[max_positive_peak_idx]/region_signal[first_negative_peak_idx])
        
        rising_rate =(region_signal[max_positive_peak_idx]-region_signal[0]) / (max_positive_peak_idx / fs)
        
        
        maxtofirstratio=abs(region_signal[max_positive_peak_idx]/region_signal[0])
        
        
        
    elif max_positive_peak_idx is not None:
        
        ffrising_rate=None
        ffmaxtofirstratio=None
        # Rising rate between baseline (assumed 0) and first positive peak
        rising_rate = (region_signal[max_positive_peak_idx]-region_signal[0]) / (max_positive_peak_idx / fs)
        maxtofirstratio=abs(region_signal[max_positive_peak_idx]/region_signal[0])
    else:
        rising_rate = None
        maxtofirstratio=None#
        ffrising_rate=None
        ffmaxtofirstratio=None#No valid rising rate if no positive peak exists

    # Calculate falling rate
    if max_positive_peak_idx is not None and second_negative_peak_idx is not None and aftermaxbaseline is not None:
        ffalling_rate = (
            region_signal[max_positive_peak_idx] - region_signal[second_negative_peak_idx]
        ) / ((second_negative_peak_idx - max_positive_peak_idx) / fs)
        
        fmaxtosecratio=abs(region_signal[max_positive_peak_idx]/region_signal[second_negative_peak_idx])
        
        falling_rate = (
            region_signal[max_positive_peak_idx] - aftermaxbaseline
        ) / ((aftermax_idx-max_positive_peak_idx) / fs)
        
        srisingrate= (region_signal[second_negative_peak_idx] -  aftersecbaseline
        ) / ((second_negative_peak_idx - aftersec_idx) / fs)
        
        smaxtosecratio=abs(region_signal[second_negative_peak_idx]/region_signal[aftersec_idx])
        
        
    elif max_positive_peak_idx is not None and aftermaxbaseline is not None:
        falling_rate = (
            region_signal[max_positive_peak_idx] - aftermaxbaseline
        ) / ((aftermax_idx-max_positive_peak_idx) / fs)
        
        fmaxtosecratio=None
        ffalling_rate=None
        srisingrate=None
        smaxtosecratio=None
        
        
        
    else:
        falling_rate=None
        fmaxtosecratio=None
        ffalling_rate=None
        srisingrate=None
        smaxtosecratio=None
        
   
        
        
    # No valid falling rate if peaks are missing

    # Peak-to-Peak Amplitude
    peak_to_peak_amplitude = max_positive_peak - max_negative_peak

    # Peak-to-Peak Durations
    first_to_last_negative_duration = (
        (second_negative_peak_idx - first_negative_peak_idx) / fs
        if first_negative_peak_idx is not None and second_negative_peak_idx is not None
        else None
    )
    positive_to_last_negative_duration = (
        (second_negative_peak_idx - max_positive_peak_idx) / fs
        if second_negative_peak_idx is not None and max_positive_peak_idx is not None
        else None
    )

    # Results
    analysis_result = {
        "Rising Rate": rising_rate,
        "Falling Rate": falling_rate,
        "Max Peak Rising Rate": ffrising_rate,
        "Second Peak Rising Rate": srisingrate,
        "Max Peak Falling Rate": ffalling_rate,
        "Max Positive Peak": max_positive_peak,
        "Max Positive Peak Index": max_positive_peak_idx,
        "Max Negative Peak": max_negative_peak,
        "Max Negative Peak Index": max_negative_peak_idx,
        "Peak-to-Peak Amplitude": peak_to_peak_amplitude,
        "First to Last Negative Peak Duration (s)": first_to_last_negative_duration,
        "Positive to Last Negative Peak Duration (s)": positive_to_last_negative_duration,
        "First Positive Peak": region_signal[first_positive_peak_idx] if first_positive_peak_idx is not None else None,
        "First Positive Peak Index": first_positive_peak_idx,
        "First Negative Peak": region_signal[first_negative_peak_idx] if first_negative_peak_idx is not None else None,
        "First Negative Peak Index": first_negative_peak_idx,
        "Second Negative Peak": region_signal[second_negative_peak_idx] if second_negative_peak_idx is not None else None,
        "Second Negative Peak Index": second_negative_peak_idx,
        "Max To baseline peak ratio":  maxtofirstratio,
        "Second To Baseline peak ratio": smaxtosecratio, 
        "Max To first peak ratio":  ffmaxtofirstratio,
        'Trace':[trace]
        
    }

    # Plot the region with the selected peaks

    if plotting==True:

        x_values = np.arange(combined_region["Start Index"], combined_region["End Index"])
        plt.figure(figsize=(12, 6))
        plt.plot(x_values, region_signal, label="Region Signal", color="blue")
        plt.axhline(y=threshold)
        plt.axhline(y=-1*threshold)
        
        if  max_positive_peak_idx !=None:
            plt.scatter(x_values[max_positive_peak_idx], region_signal[max_positive_peak_idx], 
                        color="green", label="Max Positive Peak")
            
        if  first_negative_peak_idx !=None:
            plt.scatter(x_values[first_negative_peak_idx], region_signal[first_negative_peak_idx], 
                        color="red", label="First Negative Peak")
            
        if  second_negative_peak_idx !=None:
        
            plt.scatter(x_values[second_negative_peak_idx], region_signal[second_negative_peak_idx], 
                        color="brown", label="Second Negative Peak")
        plt.title("Selected Region with Extreme Peaks")
        plt.xlabel("Sample Index")
        plt.ylabel("Amplitude")
        plt.legend()
        plt.grid(True)
        plt.show()

    return analysis_result


# In[ ]:




