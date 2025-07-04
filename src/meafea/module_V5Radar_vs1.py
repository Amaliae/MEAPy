#!/usr/bin/env python
# coding: utf-8

# In[1]:


from statistics import mean
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy.stats import wilcoxon
from scipy import stats

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

import neo

import quantities as pq
from elephant import statistics
from elephant import kernels

from scipy.stats import wilcoxon


# In[6]:


def quality(dataset, Param):
    
    
    ##normal or not if not what else
    #try modification , whether it helps
    #shoulbe filtereed
    
    dataset=Param
    
    


# In[7]:


def load_raw(path, identifier, identifier2):
    
    dirs = os.listdir(path)
    
    
    filename=[]
    
   
    
    
    
    for f in dirs:
        
        if re.search('h5', f) and re.search(identifier, f) and re.search(identifier2, f):
            
            filename.append(f)
            
    
            
            
    return filename


# In[8]:


def led_info(Recordingfilename, Metadatafile):
    
    """Fuction to extract the LED stimulated well label, led color, led times and duration of the pulse"""
    
    R=McsPy.McsData.RawData(Recordingfilename)
    
    
    Info={}
    
   
    
    
    
    LED_T=[]
    LED_D=[]
    
    for i in list(R.recordings[0].event_streams[3].event_entity.keys()):
        
        Wells=R.recordings[0].event_streams[3].event_entity[i].info.__dict__['info']['SourceChannelLabels']
        Label=R.recordings[0].event_streams[3].event_entity[i].info.__dict__['info']['Label']
        LED_T.extend(R.recordings[0].event_streams[3].event_entity[i].data[0].tolist())
    
        LED_D.extend(R.recordings[0].event_streams[3].event_entity[i].data[1].tolist())
    
    Info['Wells']=Wells
    Info['Label']=Label
    Info['Time']=LED_T
    Info['Duration']=LED_D
  
    
    return Info


# In[9]:


def load_exported0(path, identifier, identifier2):
    
    dirs = os.listdir(path)
    
    df_bursts=0
    df_spikes=0
    
    
    filename=[]
    
    i=0
    j=0
    
   
    
    for f in dirs:
        
        #print(f)
        #print(f[:-10])
       
        
        if re.search('.csv$', f) and re.search(identifier, f) and re.search(identifier2, f):
            
            if re.search('Spikes', f):
                
                if i==0:
                    
                    df_spikes=pd.read_csv(path+'/'+f)
                    df_spikes.drop(['Unnamed: 0'], axis=1, inplace=True)
                    
                    df_spikes['File']=f
                    #df_spikes['Experiment']=df_spikes['Experiment'].apply(lambda x : f[:-10])
                    #df_spikes['Dose Label']=df_spikes['Dose Label'].apply(lambda x : str(x))

                    i=i+1

                else:
                    to_c=pd.read_csv(path+'/'+f)
                    to_c['File']=f
                    to_c.drop(['Unnamed: 0'], axis=1, inplace=True)
                    #to_c['Experiment']=to_c['Experiment'].apply(lambda x : f[:-10])
                    #to_c['Dose Label']=to_c['Dose Label'].apply(lambda x: str(x))
                    df_spikes=pd.concat([df_spikes, to_c])
                    
                
               
            
            if re.search('Bursts', f):
                
                
                
                if j==0:
                    
                    df_bursts=pd.read_csv(path+'/'+f, index_col=[0])
                    df_bursts['File']=f
                    #df_bursts['Experiment']=df_bursts['Experiment'].apply(lambda x : f[:-10])  
                    j=j+1
                
                else:
                    to_c=pd.read_csv(path+'/'+f, index_col=[0])
                    to_c['File']=f
                    #to_c['Experiment']=to_c['Experiment'].apply(lambda x : f[:-10]) 
                    
                    df_bursts=pd.concat([df_bursts,to_c])
                    
                
                
            filename.append(f)
            
          
    
    
   
    
  
    try:
        A=dict(zip(np.unique(df_spikes['Channel Label'].values), np.arange(0, 12, 1)))
        df_spikes['True Label']=df_spikes['Channel Label']
        df_spikes['Channel Label']=df_spikes['Channel Label'].apply(lambda x : A[x])
        
        
        
    except: 
        'No spike files'
    try:
        A=dict(zip(np.unique(df_bursts['Channel Label'].values), np.arange(0, 12, 1)))
        df_bursts['True Label']=df_bursts['Channel Label']
        df_bursts['Channel Label']=df_bursts['Channel Label'].apply(lambda x : A[x])
        
        
    except:
        'No burst file'
    try:
        
        df_spikes['Dose Label']=df_spikes['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)
    except: 
        'No spike files'
    
    try:
        df_bursts['Dose Label']=df_bursts['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)

        columnsburst=df_bursts.columns

        columnsburst=['Timestamp [µs]' if column=='Start timestamp [µs]' else column for column in columnsburst]

        #print(columnsburst, 'columnsb')

        df_bursts.columns=columnsburst
    except:
        'No burst file'
        
        
    
    
    
        
    return (df_spikes, df_bursts, filename)


# In[32]:


def load_exported(path, identifier, identifier2, newpharm, Comps):
    
    
    
    dirs = os.listdir(path)
    
    df_bursts=0
    df_spikes=0
    
    
    filename=[]
    
    i=0
    j=0
    
   
    
    for f in dirs:
        
        #print(f)
        #print(f[:-10])
       
        
        if re.search('.csv$', f) and re.search(identifier, f) and re.search(identifier2, f):
            
            if re.search('Spikes', f):
                
                if i==0:
                    
                    df_spikes=pd.read_csv(path+'/'+f)
                    df_spikes.drop(['Unnamed: 0'], axis=1, inplace=True)
                    
                    df_spikes['File']=f
                    #df_spikes['Experiment']=df_spikes['Experiment'].apply(lambda x : f[:-10])
                    #df_spikes['Dose Label']=df_spikes['Dose Label'].apply(lambda x : str(x))

                    i=i+1

                else:
                    to_c=pd.read_csv(path+'/'+f)
                    to_c['File']=f
                    to_c.drop(['Unnamed: 0'], axis=1, inplace=True)
                    #to_c['Experiment']=to_c['Experiment'].apply(lambda x : f[:-10])
                    #to_c['Dose Label']=to_c['Dose Label'].apply(lambda x: str(x))
                    df_spikes=pd.concat([df_spikes, to_c])
                    
                
               
            
            if re.search('Bursts', f):
                
                
                
                if j==0:
                    
                    df_bursts=pd.read_csv(path+'/'+f, index_col=[0])
                    df_bursts['File']=f
                    #df_bursts['Experiment']=df_bursts['Experiment'].apply(lambda x : f[:-10])  
                    j=j+1
                
                else:
                    to_c=pd.read_csv(path+'/'+f, index_col=[0])
                    to_c['File']=f
                    #to_c['Experiment']=to_c['Experiment'].apply(lambda x : f[:-10]) 
                    
                    df_bursts=pd.concat([df_bursts,to_c])
                    
                
                
            filename.append(f)
            
          
    
    
   
    
  
    try:
        A=dict(zip(np.unique(df_spikes['Channel Label'].values), np.arange(0, 12, 1)))
        df_spikes['True Label']=df_spikes['Channel Label']
        df_spikes['Channel Label']=df_spikes['Channel Label'].apply(lambda x : A[x])
        
        
        
    except: 
        print('No spike files')
    try:
        A=dict(zip(np.unique(df_bursts['Channel Label'].values), np.arange(0, 12, 1)))
        df_bursts['True Label']=df_bursts['Channel Label']
        df_bursts['Channel Label']=df_bursts['Channel Label'].apply(lambda x : A[x])
        
        
    except:
        print('No burst file', 'can not assign chan')
    try:
        
        df_spikes['Dose Label']=df_spikes['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)
        
        
        if newpharm==True:
            df_spikes['Dose Label']=df_spikes['Compound ID']
            
            
            df_spikes['Compound ID']=df_spikes['Well ID'].apply(lambda x : Comps[int(x)])
            
            #print(np.unique(df_spikes['Compound ID'].values), 'inloadexpoetedcompaounds')
            
            #print(np.unique(df_spikes['Compound ID'].values), 'after')
            
    except: 
        print('No spike files')
    
    try:
        
        df_bursts['Dose Label']=df_bursts['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)
        
        
        if newpharm==True:
            
            #print(df_bursts['Compound ID'].unique(), 'compound')
            
            df_bursts['Dose Label']=df_bursts['Compound ID']
            
            #print(df_bursts['Dose Label'].unique(), 'dose')
            
            

            df_bursts['Compound ID']=df_bursts['Well ID'].apply(lambda x : Comps[int(x)])
            
            
        #df_bursts['Dose Label']=df_bursts['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)

        columnsburst=df_bursts.columns
        
        #print(columnsburst, 'cols')

        columnsburst=['Timestamp [µs]' if column=='Start timestamp [µs]' else column for column in columnsburst]

        #print(columnsburst, 'columnsb')

        df_bursts.columns=columnsburst
        
        
    except:
        print('No burst file, newpharm?', Comps)
        
    
        
    
    
            
    
        
    return (df_spikes, df_bursts, filename)


# In[33]:


def experiment_ledinfo(folder_path, save_to, identifier1, identifier2, MetadataMaskPath,  dead_time=3, threshold=4.8, negative=True, positive=True, fs=20000):
    
    
    filenames=load_raw(folder_path, identifier1, identifier2)
    
    Entry={}
    
    i=0
    
    for file in filenames:
        
        i=i+1
        
        #print(file)
        
        Metadatafile=folder_path+'/'+file[:-7]+'.mws'
        
        Recordingfile=folder_path+'/'+file
        
        tree = etree.parse(Metadatafile)
        ExperimentID=tree.xpath('//ExperimentID/text()')[0]
        Compound=tree.xpath('//CompoundName/text()')[0] 
        L=led_info(Recordingfile, Metadatafile)
        Entry[ExperimentID+Compound+str(i)]=L
        
     
          
    return  Entry


# In[ ]:





# In[ ]:





# In[ ]:





# In[34]:


def bin(spikes, start, end, winsize):
    
    spikes=spikes[(spikes>start) & (spikes<end)]
    
   
    
    #winsizeopt=est_winsize(spikes)
    ar=np.zeros(int((end-0)/winsize)+2)
    #aropt=np.zeros(int((end-0)/winsizeopt)+2)
    ##print('winsizeopt', winsizeopt)
    
    
    
    
    
   
    
    for i in spikes: 
        
        

        idx=int(round(i//winsize)) 
        #idxopt=int(round(i//winsizeopt)) 
        
        
        
        
        ar[idx]=ar[idx]+1
        
        #aropt[idxopt]=aropt[idxopt]+1
        
        
    #plt.plot(ar[int(round(start/winsize)):])
    #plt.xlabel('bin N')
    #plt.ylabel('N of spikes')
    #plt.title('Bining')
    #plt.show()
    
    
    return ar[int(round(start/winsize)):]


# In[35]:


def moving_Average(bins, to, channel, ID, plot):
    
    
    
    kernel_size = len(bins)//10
    kernel = np.ones(kernel_size) / kernel_size
    bin_conv = np.convolve(bins, kernel, mode='valid')
    if plot=='True':
        plt.figure(figsize=(10, 8))
        plt.plot(bins, "-b", label="bins")
        plt.plot(bin_conv, "-r", label="bins conv")
        plt.xlabel('Bin')
        plt.ylabel('N of Spikes')
        plt.legend(loc="upper left")
        plt.savefig(to+'/'+'Rollingaverage'+str(channel)+ID+'.png')
        plt.show()
    
    return bin_conv


# In[36]:


def cumsmooth(spikes, interval, delta, start, end):
    
    time=np.zeros(end-start)
    
    time[spikes-start]=1
   
    
    cumaxis=np.cumsum(time)
    
    upper_index=np.arange(interval, len(time), delta)
    
    ##print(upper_index, 'upper_boundary')
    
    lower_index=np.arange(0, len(time)-interval, delta)
    
    ##print(lower_index,  'lower_boundary')
    
    values_upper_index=cumaxis[upper_index]
    value_lower_index=cumaxis[lower_index]
    
    smoothed=values_upper_index-value_lower_index
    
    ##print(smoothed, 'smoothed')
    
    #plt.plot(smoothed)
    
    #plt.xlabel('Interval order')
    #plt.ylabel('N of Spikes')
    #plt.title('Cumsmooth')
    
    #plt.show()
    
    
    return smoothed


# In[ ]:





# In[37]:


def mfr_bfr_label(to, ID, compound, series, lag, winsize, start, end, s, e, 
                  plotsi, label_index, ch_max, label, df_ISI, Type, acf, group_iN,  scale, ISI_group):
    
    
    
   
    
    s=int(np.zeros(int((end-0)/winsize)+2)[int(round(start/winsize)):].shape[0])
    sm=moving_Average(np.ones(int((end-0)/winsize)+2)[int(round(start/winsize)):], to, 'sham', ID, 'False').size
    
    #s=cumsmooth(np.zeros(int(end-0)), 10000, 5000, start, end)
    #integers _ 
    Nbin=len(np.arange(0, 20, 1))-1
    
    arr_b=np.zeros([24, 12, s])
    arr_isi=np.zeros([24, 12, Nbin])
    acf_b=np.zeros([24, 13, sm])
    arr_m_b=np.zeros([24, 12])
    arr_m_ch=np.zeros([24, 12])
    arr_std_b=np.zeros([24, 12])
    duration=(end-start)/1000000
    
    arr_w_b=np.zeros([24, 12])
    
    
    rasters=np.empty([12, ch_max])
    wellssttcpre=np.zeros([24, 12])
    
    #print(sm, 'sm')
   
   
    
    dict_ISI={'ISI':[], 'Group':[]}
    
    isic=[] 
    gp=[]
        
    wls=[]
    lightwave=[]
    
    #for each well
        
    
    figure_collect=[]
    
    series['Well ID']=series['Well ID'].apply(lambda x : int(x))
    wells=np.unique(series['Well ID'].sort_values(ascending=True).values)
    
    xls=start+1
    xle=end
    
    for indexwell, i in enumerate(wells):
        
        channels=np.unique(series.loc[series['Well ID']==i]['Channel Label'].sort_values(ascending=True).values)
        channelids=np.unique(series.loc[series['Well ID']==i]['Channel ID'].sort_values(ascending=True).values)
        all_spikes=len(series.loc[series['Well ID']==i])
        
        spikesynpre=[]
        labrasters=[]
        
       
        
        if len(plotsi)>0:
            
            if plotsi[indexwell][1].ndim==2:
            
        
                rasteraxes=plotsi[indexwell][1][0, label_index]

                MAaxes=plotsi[indexwell][1][1, label_index]

                logisiaxes=plotsi[indexwell][1][2, label_index]

            else:
                rasteraxes=plotsi[indexwell][1][0]

                MAaxes=plotsi[indexwell][1][1]

                logisiaxes=plotsi[indexwell][1][2]
            
        
        
        
        
        
        
        colors = ['C{}'.format(i) for i in range(0, len(channels))]
        
        colors=colors+['black']
        
        lineoffset=np.arange(0, 12, 1).tolist()
        linelenght=np.ones(12)/2
        rasterlab=[str(i) for i in range(12)]

       
        
        for ch in range(len(channels)):
            
            
            
            if Type=='Spike':
                
                print(i, 'well')
            
                spikes_c=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label) & 
                        (series['Well ID']==i)]['Timestamp [µs]'].sort_values(ascending=True).values.astype('int64')
                ChiD=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label) 
                                & (series['Well ID']==i)]['Channel ID'].values[0]
                
                print(ChiD, 'Channel ID')
                
            elif Type=='Burst':
                
                spikes_c=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)]['Timestamp [µs]'].sort_values(ascending=True).values.astype('int64')
                ChiD=series.loc[(series['Channel Label']==channels[ch])  & (series['Well ID']==i)]['Channel ID'].values[0]

            b=bin(spikes_c, start, end, winsize)
            
            #g=cumsmooth(spikes_c, 2000000, 500000, start, end)
            #spiketrain = neo.SpikeTrain(spikes_c/1000, t_start=(start/1000)*pq.ms, t_stop=(end/1000)*pq.ms, units='ms')

            #kernel = kernels.GaussianKernel(sigma=1000 * pq.ms)

            #rate = statistics.instantaneous_rate(spiketrain, sampling_period=20*pq.ms, kernel=kernel)
            
           
            
            ##print(rate.shape, ch, 'rate shape')
            
            spikes_c=spikes_c[(spikes_c>start) & (spikes_c<end)]
            
            
            
            spikesynpre.append(spikes_c)
        
    
            
            
            
            bins_conv=moving_Average(b, to, str(ch)+" "+"W"+str(i), ID, 'False')
            
            #bins_conv=slowess.nonparametric.lowess(exog=np.arange(0, len(b), 1), endog=b, frac=0.3)[:, 1]
            
            #coefs=np.polyfit(np.arange(0, len(b), 1), b, 2)
            #bins_conv=np.polyval(coefs, b)
            
        
            acf_b[int(i), int(ch), :]=bins_conv
            
            
            
            
            
            #bin
            
            #Raster plor per channel
            
            to_pad=np.zeros(int(ch_max-len(spikes_c))).tolist()  #padding for rasters
            spike=np.array(to_pad+spikes_c.tolist())
            rasters[channels[ch], :]=spike
           

            arr_b[int(i), int(ch), :]=b
            labrasters.append(ch)
            #acf_b[int(i), int(ch), :]=scipy.signal.correlate(b, b)[:lag] 
            
           

           #log 10 for firing rate per channel 
            arr_m_b[int(i), int(ch)]=logMFR(spikes_c, duration)  #/((all_spikes)/duration))*logMFR(spikes_c, duration)#np.mean(b)*(1000000/(winsize))
            arr_m_ch[int(i), int(ch)]=ChiD
           
            
            #standart deviation for mean firing rate
            arr_std_b[int(i), int(ch)]=np.std(b)*(1000000/(winsize)) 
            
            #interspike intervals in a single channel 
            isc=np.log10(np.diff(np.array(spikes_c)))
            
            isicount, isibins=np.histogram(isc, bins=np.arange(0, 20, 1))
            
            isibins=[i+np.diff(isibins)[0] for i in isibins][:-1]
            
            arr_isi[int(i), int(ch), :]=isicount
            
            if len(plotsi)>0:
                
            
                logisiaxes.plot(isibins, isicount, color=colors[ch])
                logisiaxes.set_xlabel('Bin center', fontsize=55)
                logisiaxes.set_ylabel('Count', fontsize=55)
                logisiaxes.tick_params(axis='x', labelsize=40) 
                logisiaxes.tick_params(axis='y', labelsize=40)

                MAaxes.plot(bins_conv, color=colors[ch])

                MAaxes.set_xlabel('Bins', fontsize=55)
                MAaxes.set_ylabel('N of Spikes', fontsize=55)


                MAaxes.tick_params(axis='x', labelsize=40) 
                MAaxes.tick_params(axis='y', labelsize=40)

            
            isc=isc[isc>=0]   #ISi

            isic.extend(isc) #extend for each channel 

            gp.extend(np.repeat(int(i), len(isc)).tolist())  #group name extension 
            wls.extend(np.repeat(group_iN[int(i)], len(isc)).tolist())
            lightwave.extend(np.repeat('black', 12).tolist())
            
        if len(plotsi)>0:
            rasteraxes.eventplot(rasters, lineoffsets=lineoffset, colors='black', linelengths=linelenght, linewidths=0.3)
            rasteraxes.set_xlim(xls, xle)
            #rasteraxes.set_yticklabels(np.arange(rasters.shape[0]), rasterlab, fontsize='x-small')
            rasteraxes.set_xlabel('Time', fontsize=55)
            rasteraxes.set_title(label, fontsize=70)

            rasteraxes.set_ylabel('Channel', fontsize=55)
            rasteraxes.tick_params(axis='x', labelsize=40) 
            rasteraxes.tick_params(axis='y', labelsize=40)
            
            rasteraxes.set_xticks(np.linspace(xls, xle, 5))
            
            rasteraxes.set_xticklabels(np.linspace(xls//1000000, xle//1000000, 5))

            MAaxes.plot(np.mean(acf_b[int(i), :, :], axis=0), color=colors[-1])
            logisiaxes.plot(isibins, np.mean(arr_isi[int(i), :, :], axis=0), color=colors[-1])
            
        if len(plotsi)>0:
            
            figrast, axrast=plt.subplots(1, 1, figsize=(16, 10))

            axrast.eventplot(rasters, lineoffsets=lineoffset, colors='black', linelengths=linelenght, linewidths=0.5)
            axrast.set_xlim(xls+11000000, xls+41000000)
            axrast.set_yticklabels(labrasters)
            axrast.set_xlabel('Time', fontsize=55)
            axrast.set_title(label,  fontsize=50)
            axrast.set_xticks(np.linspace(xls+11000000, xls+41000000, 5))

            axrast.set_xticklabels(np.linspace((xls+11000000)//1000000,  (xls+41000000)//1000000, 5))


            axrast.set_ylabel('Channel', fontsize=55)
            axrast.tick_params(axis='x', labelsize=40) 
            axrast.tick_params(axis='y', labelsize=40)

            #figrast.savefig(to+'/'+'rasters'+label+compound+ID+str(i)+'.png')
        
        
        
        presttcs=Channel_STTC(spikesynpre, start, end)
       
        
        wellssttcpre[int(i), :]=presttcs 
    


        
            #extend for each group
        if acf==True:

            plt.figure()
        
            

            mean_b=np.mean(arr_b[int(i), :, :], axis=0)

            plot_acf(mean_b, lags=np.arange(1, 500, 1), alpha=0.05)

            plt.title('AUTOCORR')
            plt.savefig(to+'/'+str(i)+label+ID+'ACF'+Type+str(i)+str(ch)+'.tif', format='tif', bbox_inches='tight')
            
            plt.figure(figsize=(10, 8))
            plt.plot(np.mean(acf_b[int(i), :, :], axis=0), "-b", label="bins")
            
            plt.xlabel('Bin')
            plt.ylabel('N of Spikes')
            plt.legend(loc="upper left")
            plt.savefig(to+'/'+'Rollingaverage'+str(i)+ID+'.png')
            plt.show()
            
        
          
    #ISI append metadata 
    dict_ISI['Well ID']=gp
    dict_ISI['ISI']=isic
    dict_ISI['Group']=wls
    I=pd.DataFrame(dict_ISI)
    
    I['Stim']=np.repeat('off', len(I))
    I['Experiment']=np.repeat(ID, len(I))
    I['Compound']=np.repeat(compound, len(I))
    I['Label']=np.repeat(label, len(I))
    
    I['StimWells']=np.repeat('off', len(I))
    
    plt.figure(figsize=(10, 8))
    plt.plot(np.mean(acf_b, axis=(0, 1)), "-b", label="bins")

    plt.xlabel('Bin')
    plt.ylabel('N of Spikes')
    plt.legend(loc="upper left")
    plt.savefig(to+'/'+'Rollingaverage'+label+compound+ID+'.png')
    plt.show()
    
    #figrast, axrast=plt.subplots(1, 1, figsize=(16, 10))
    
    #axrast.plot(np.linspace(0, 20, acf_b.shape[2]), np.mean(np.mean(acf_b, axis=1), axis=0))

     
    #figrast.savefig(to+'/'+'Binned Average of plate'+label+compound+ID+'.png')
        

    
    
    
   
    
    
    
    
    return (np.mean(arr_isi, axis=1), arr_m_b, I, wells,  arr_m_ch, lightwave, plotsi, [], [], [], [], [], [wellssttcpre], [])


            


# In[38]:


#function correlated wells 

#different channels in a single well are in the different phases, or correlation 

#if not correlations, 
# determine phase of the periodic signal (each well in different phase)
#and for the each group peak and compare the bins, in the same phase (angle)




#plan fft transform of the data: 
#summed power of the given range, and then comparison of the ranges: 


#so then, dataframe results and concat


# In[ ]:





# In[39]:


def FFT(conved_bin, winsize, duration, ranges):
    from scipy.fft import rfft, rfftfreq
    
    
    results={}
    
    
    #normalized conved_bin, this is  a numpy array
    
    signal=conved_bin/conved_bin.max()
    
    #winsize, duration in us
    
    samplingrate=1/winsize(1000000)
    
    duration=duration/1000000
    
    N = samplingrate * duration
    
    yf = rfft(signal)  #power
    xf = rfftfreq(N, 1 / samplingrate) #frequencies
    
    #ranges are user deined list of tuples, with min and max values of ranges 
    for i in ranges:
        
        if i[1]<samplingrate/2:
            
            indexes=np.where((xf>i[0])and (xf<i[1]))[0]
            
            amplitudes=np.sum(np.abs(yf[indexes]))
            
            results[i]=amplitudes
            
            
    
    return results
    
    
    
    
    


# In[40]:


def fake_led(times, start, end):
    
    diff=np.diff(times)[0]
    
    se=np.arange(times[0], start, -1*diff)[::-1][:-1]
    ee=np.arange(times[-1], end, diff)[1:]
    
    
    return(se.tolist(), ee.tolist())
       


# In[41]:


def ledprepost(chspikes, st, et, times):
    
    fake=[]
    realpre=[]
    realpost=[]
    
    chspikes=np.array(sorted(chspikes))
    
    ###those already spikes outside of tange 
    
    for t in (st+times+et): 
        
        spkpre=chspikes[chspikes<t]
        
        spkpos=chspikes[chspikes>t]
        
        
        spkpret=t-spkpre
        
        spkpost=spkpos-t
        
        realpre.append(spkpre[-10:])
        
        realpost.append(spkpos[:10])
        
        
    return (realpre, realpost)


# In[42]:


def diff_prepostled(spikecounts, st ,et, times):
    
    preled=spikecounts[0]
    postled=spikecounts[1]
    
    ##Idon't need each to have 10 preposit spikes, but i need them to be equal
    
    ##padding with range limit values 
    
    #preledpadded=[np.ones([10-len(preled[l])])*prer.tolist()+preled[l]  for l in range(len(preled))]
    
    #postledpadded=[postled[l]+np.ones([10-len(postled[l])])*postr.tolist()  for l in range(len(postled))]
    
    
    #preledc=np.array([preledpadded])

    #postledc=np.array([postledpadded])
    
    
    preledind=[l  for l in range(len(preled)) if len(preled[l])==10]
    
    
    ##10 spikes stimuli


    postledind=[l for l in range(len(postled)) if len(postled[l])==10]

    commonind=np.intersect1d(preledind, postledind)
    
    ##3take only those stimuli

    preledc=np.array([preled[l] for l in commonind])

    postledc=np.array([postled[l] for l in commonind])
    
    
    if (preledc.shape[0])>0: 
        
        
        allleds=np.concatenate([preledc, postledc], axis=1)
        
        fakes=np.array(st+times.tolist()+et)

        fakesin=fakes[commonind]

        commonintrue=[cm for cm in range(len(commonind)) if np.isin(fakes[commonind[cm]], times)==True]

        commoninfalse=[cm for cm in range(len(commonind)) if np.isin(fakes[commonind[cm]], times)==False][:len(commonintrue)]

        allledstrue=np.diff(allleds[commonintrue, :]/1000, axis=1)

        allledsFalse=np.diff(allleds[commoninfalse, :]/1000, axis=1)
        
        low25fake=np.percentile(allledsFalse[:, 9], 25)
        
        succstim=np.where(allledstrue[:, 9]<low25fake)[0]
        
        print(allledstrue.shape[0])
        
        
        return (allledstrue, allledsFalse, succstim)
        
        
        


# In[43]:


def FFT(conved_bin, winsize, duration, ranges):
    from scipy.fft import rfft, rfftfreq
    
    
    results={}
    
    
    #normalized conved_bin, this is  a numpy array
    
    signal=conved_bin/conved_bin.max()
    
    #winsize, duration in us
    
    samplingrate=1/winsize(1000000)
    
    duration=duration/1000000
    
    N = samplingrate * duration
    
    yf = rfft(signal)  #power
    xf = rfftfreq(N, 1 / samplingrate) #frequenciesnp_fft = np.fft.fft(ys)
    #amplitudes = 2 / n_samples * np.abs(np_fft) 
    #frequencies = np.fft.fftfreq(n_samples) * n_samples * 1 / (t1 - t0)
    
    
    #np_fft = np.fft.fft(ys)
    #amplitudes = 2 / n_samples * np.abs(np_fft) 
    #frequencies = np.fft.fftfreq(n_samples) * n_samples * 1 / (t1 - t0)
    
    #ranges are user deined list of tuples, with min and max values of ranges 
    for i in ranges:
        
        if i[1]<samplingrate/2:
            
            indexes=np.where((xf>i[0])and (xf<i[1]))[0]
            
            amplitudes=np.sum(np.abs(yf[indexes]))
            
            results[i]=amplitudes
            
            
    
    return results
    
    
    


# In[44]:


def find_amp_weird(ID, listamps):
    
    amplitude=14
    
    try:
    
        amp_idx=ID.find('A')

        for amp in listamps:
            if str(amp) in ID[amp_idx-2:amp_idx+3]:

                amplitude=amp
                
    except:
        'no info'
            
    return amplitude
    
   
    
   
        


# In[ ]:





# In[ ]:





# In[1]:


def label_to_id_old(led_label, WellID, Welllabel):
    
    #the led and stimulation are given by well labels, we work on ids. 
    led_id=[]
    
    led_label=led_label.split(',')
    l_condition=np.repeat('off', 24)
    for label in led_label:
        
        
        
        i=int(WellID[np.where(Welllabel==label)[0][0]])
        l_condition[i]='On'
        led_id.append(i)
        
    return led_id, l_condition
        


# In[46]:


#mfr_bfr_Stim (LED_dict.key() contain Label)


# In[47]:


def mfr_bfr_Stim(to, ID, compound, series, lag, winsize, start, end, s, e, 
                  plotsi, label_index, ch_max, label, df_ISI, Type, Param, acf, group_iN, LED_dict, ledid,  ISI_group, osuc,
                 listamp):
    
   
    
    
    s=int(np.zeros(int((end-0)/winsize)+2)[int(round(start/winsize)):].shape[0])
    
    #print(s, 's')
    
    kernel_shape=moving_Average(np.ones(int((end-0)/winsize)+2)[int(round(start/winsize)):], to, 'sham', ID, 'False').size
    
    
    #shamtrain=np.arange(start/1000, end/1000, 0.05)
    
   
    #shamtrainneo=neo.SpikeTrain(shamtrain, t_start=(start/1000)*pq.ms, t_stop=(end/1000)*pq.ms, units='ms')

    #kernel=kernels.GaussianKernel(sigma=2000 * pq.ms)

    #rate = statistics.instantaneous_rate(shamtrainneo, sampling_period=0.05*pq.ms, kernel=kernel)
    
    #sm=rate.reshape(rate.shape[0]).size
    
   
    #moving_Average(rate.reshape(rate.shape[0]), to, 'sham', ID, 'False').size
   
 
     
    
    ##Averaging instanteous FR
    arr_m_b=np.zeros([24, 12])
    arr_m_ch=np.zeros([24, 12]) 
    arr_b=np.zeros([24, 13, s])
    light_corr=np.zeros([24, 12])
    light_lag=np.zeros([24, 12])
    
    wellssttcpre=np.zeros([24, 12])
    wellssttcpost=np.zeros([24, 12])
    wellssttclight=np.zeros([24, 12])
    
   
    
    
    #placeholders for all wells and electrodes
    
    ###to extract Channel IDs
    
    duration=(end-start)/1000000  ###duration
    
    
    arr_w_b=np.zeros([24, 12])  ###? 
    
    
    rasters=np.empty([13, ch_max])
    ###rasters in Stim function always has LED axis
    
   
    
    dict_ISI={'ISI':[], 'Group':[]}
    
    isic=[] 
    gp=[]
        
    wls=[]
    
    light=[]
    I={}
    lightwave=[]
    #for each well
        
  
    
    series['Well ID']=series['Well ID'].apply(lambda x : int(x))
    wells=np.unique(series['Well ID'].sort_values(ascending=True).values)
    #print('wells' , wells)
    
    xls=start+1
    xle=end
    
    Rasters=np.empty([len(wells)+1, ch_max])
    
    
    ####LED####
    
       
    WellID=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellID.csv.npy")
    Welllabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellLabel.csv.npy")
    Channellabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\ChannelLabel.csv.npy")
    
    
    if 'Label' in LED_dict.keys():
        
        LED_IDs=[]
        
        
        
        # this is LED metadata 
        
        for lw in LED_dict['Wells']:
            
        
            LED_labels=lw
            
            #extracts and modifies stimulated wells, also modify this
            led_wells=label_to_id(LED_labels, WellID, Welllabel)[0]
            
            LED_IDs.extend(led_wells)
        
        light_label=LED_dict['Label']
        
        light_label=[string[-6:].lower() for string in light_label]
        
        #print(light_label, 'light_label')
        

        
       
        
        light_dur=LED_dict['Duration'][0][0]/1000


        LED=np.unique(np.array(LED_dict['Time'][0]))
        
        light_freq=(1000000*len(LED))/(LED[-1]-LED[0])
        
        
        if len(listamp)>0:
            
        
        
            amps=find_amp_weird(ID, listamp)
        else:
            amps=14
            
            #in the LED stim there is no amp 
          
        
    else:   ###electrical stimulation metadata 
        
        LED=LED_dict['Times'+str(label)]*50 #electrical stimulation 
        
        #print("LED Times", LED)# This is the stimulation already in times#change to [str(label)]
        LED_IDs=LED_dict['Wells'] #well id should be given
        LED_IDs=[int(i) for i in LED_IDs]
        amps=np.array([int(i) for i in LED_dict['Amplitudes']])   #.split(' ')])
        amps=amps[amps<0] #negative phase of biphasic stimulus
    if LED[(LED>=start) & (LED<=end)].shape[0]>0:   #the stimulation series for the given dose
        
        LED_s=int(LED[(LED>=start) & (LED<=end)][0])
        
        LED_e=int(LED[(LED>=start) & (LED<=end)][-1])
        
        LED_a=LED[(LED>=start) & (LED<=end)]
        
    else: 
        
        LED_s=end
        LED_e=end
        LED_a=[]
        
    
    s_pre=int(np.zeros(int((LED_s-0)/winsize)+2)[int(round(start/winsize)):].shape[0])
    s_led=int(np.zeros(int((LED_e-0)/winsize)+2)[int(round(LED_s/winsize)):].shape[0])
    s_post=int(np.zeros(int((end-0)/winsize)+2)[int(round(LED_e/winsize)):].shape[0])
    
 
  
    try:
        
        led_axis=np.array(np.zeros(int(ch_max-len(LED_a))).tolist()+LED_a.tolist()) 
        # construction of the stimulation axis
        
    except:  ###ch_max is the latest spike in whole plate, possible that LED is latest
        
        'Dimension error'
        ch_max=len(LED_a) #  otherwise the shape of stim as the shapes for the raster should be similiar. 
        led_axis=np.array(np.zeros(int(ch_max-len(LED_a))).tolist()+LED_a.tolist())
        rasters=np.empty([13, ch_max])
        Rasters=np.empty([len(wells)+1, ch_max])
   
     
    d2=(LED_e-LED_s)/1000000  #duration of the period in seconds 
    d1=(LED_s-start)/1000000
    d3=(end-LED_e)/1000000
    #print(d3, 'duration')
    
    arr_m_pre=np.zeros([24, 12])
    arr_m_led=np.zeros([24, 12])
    arr_m_post=np.zeros([24, 12])
    
    colorsR=[]
     
    labRasters=[]
    lineoffset1=np.arange(0, 25, 1).tolist()
    linelenght1=np.ones(25)/2
    rasterlab1=[str(i) for i in range(25)]
    
    sizepsth0=int(np.diff(LED_a).min()//5000)+int(50000//5000)+1
    sizepsth=len(gen_poin(np.diff(LED_a).min()/1000, 10))
    allpsth=np.zeros([24, 12, sizepsth0])
    

    
    
    for indexwell, i in enumerate(wells):
        
        labrasters=[]
        labRasters.append(str(i))
        
        lineoffset=np.arange(0, 13, 1).tolist()
        linelenght=np.ones(13)/2
        rasterlab=[str(i) for i in range(13)]
        
        channels=np.unique(series.loc[series['Well ID']==i]['Channel Label'].sort_values(ascending=True).values)
        channelids=np.unique(series.loc[series['Well ID']==i]['Channel ID'].sort_values(ascending=True).values)
        all_spikes=len(series.loc[series['Well ID']==i])
        
        colors = ['C{}'.format(i) for i in range(0, len(channels))]
        
        colors=colors+['black']
        
        
        
        #figpsth=plt.figure(figsize=(10, 8))
        #axpsth=figpsth.add_subplot(111)
        
        if len(plotsi)>0:
            
            if plotsi[indexwell][1].ndim==2:
            
        
                rasteraxes=plotsi[indexwell][1][0, label_index]

                MAaxes=plotsi[indexwell][1][1, label_index]

                logisiaxes=plotsi[indexwell][1][2, label_index]
                
                axpsth=plotsi[indexwell][1][3, label_index]
                axpsthlast90=plotsi[indexwell][1][4, label_index]

            else:
                rasteraxes=plotsi[indexwell][1][0]

                MAaxes=plotsi[indexwell][1][1]

                logisiaxes=plotsi[indexwell][1][2]
                
                axpsth=plotsi[indexwell][1][3]
                axpsthlast90=plotsi[indexwell][1][4]
                
                
                
                
        
        
        kernel_mfr=np.zeros(kernel_shape)
        
        wellpsth=np.zeros([12, sizepsth])
        
        spikesynpre=[]
        spikesynstim=[]
        spikesynpost=[]
        
        
        
        
        
        
        for ch in range(len(channels)):
            if Type=='Spike':
            
                spikes_c=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label) & 
                        (series['Well ID']==i)]['Timestamp [µs]'].sort_values(ascending=True).values
                ChiD=series.loc[(series['Channel Label']==channels[ch]) &  (series['Well ID']==i)]['Channel ID'].values[0]
                
            elif Type=='Burst':
                
                spikes_c=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)]['Timestamp [µs]'].sort_values(ascending=True).values
                
                ChiD=series.loc[(series['Channel Label']==channels[ch])  & (series['Well ID']==i)]['Channel ID'].values[0]

            if len(spikes_c)>0:
               
                print(len(spikes_c),  spikes_c[0], spikes_c[-1], 'spike_c', label)
            
            #print(start, end, 'start, end')
            b=bin(spikes_c, start, end, winsize)
            
           
            
            arr_b[int(i), int(channels[ch]), :]=b
            
            #spiketrain = neo.SpikeTrain(spikes_c/1000, t_start=(start/1000)*pq.ms, t_stop=(end/1000)*pq.ms, units='ms')

            #kernel=kernels.GaussianKernel(sigma=2000 * pq.ms)

            #rate = statistics.instantaneous_rate(spiketrain, sampling_period=0.05*pq.ms, kernel=kernel)
            
        
            #bins_conv=rate.reshape(rate.shape[0])
            bins_conv=moving_Average(b, to, str(ch)+" "+"W"+str(i), ID, 'False')
          
            
            #bins_conv=slowess.nonparametric.lowess(exog=np.arange(0, len(b), 1), endog=b, frac=0.3)[:, 1]
            
            #coefs=np.polyfit(np.arange(0, len(b), 1), b, 2)
            #bins_conv=np.polyval(coefs, b)
            
            
            welllabellight=Welllabel[int(ChiD)]
        
            
            kernel_mfr=np.array(kernel_mfr)+bins_conv
            
            #print(ch_max, spikes_c, int(ch_max-len(spikes_c)))
            
            to_pad=np.zeros(int(ch_max-len(spikes_c))).tolist()
            
            spike=np.array(to_pad+spikes_c.tolist())
           
            rasters[channels[ch], :]=spike
            labrasters.append(channels[ch])
           

            
            #acf_b[int(i), int(ch), :]=scipy.signal.correlate(b, b)[:lag] 
            
           

           #log 10 for firing rate per channel 
            arr_m_b[int(i), int(ch)]=(logMFR(spikes_c, duration))
            
            ######/((all_spikes)/duration))*logMFR(spikes_c, duration)#np.mean(b)*(1000000/(winsize))
            
            
            arr_m_ch[int(i), int(ch)]=ChiD
           
            
            #standart deviation for mean firing rate
            ##arr_std_b[int(i), int(ch)]=np.std(b)*(1000000/(winsize)) not needed 
            
            #interspike intervals in a single channel 
            
            
           
            
            st, et=fake_led(LED_a, start, end)
            
            
            led10=ledprepost(spikes_c, st, et, LED_a.tolist())
            
            
            ledtrue, ledfake, ledsuc=diff_prepostled(led10, st, et, LED_a)
            
            
            
            
            succLED_a=LED_a[ledsuc]
            
            
            
            
            
            
            ##ledsuc is successfull train 

           
        
            
            
            ###LED channel ####
            
            spikes_led=spikes_c[(spikes_c>LED_s) & (spikes_c<LED_e)]
            
            spikesynstim.append(spikes_led)
            
            spike_base=spikes_c[spikes_c<LED_s]
            
            spikesynpre.append(spike_base)
            
            spike_post=spikes_c[spikes_c>LED_e]
            spikesynpost.append(spike_post)

            arr_m_pre[int(i), channels[ch]]=logMFR(spike_base, d1)
            
            arr_m_led[int(i), channels[ch]]=logMFR(spikes_led, d2)#np.mean(lb)*(1000/(winsize/1000))
            arr_m_post[int(i), channels[ch]]=logMFR(spike_post, d3)#np.mean(b_post)*(1000/(winsize/1000))
            
            ledx1=int(round(LED_s/winsize))-int(round(start/winsize))-5
            ledx2=int(round(LED_e/winsize))-int(round(start/winsize))
                



            ###binned activity for this 3 periods
            b_pre=bin(spike_base, start, LED_s, winsize)  #before stimulation 
            b_post=bin(spike_post, LED_e, end, winsize)# after stimulation 
            lb=bin(spikes_led, LED_s, LED_e, winsize) #during stim
            led_b=bin(LED_a, LED_s, LED_e, winsize) #binned light
            
            synchrony=Light_corr(lb, led_b)
            light_corr[int(i), channels[ch]]=max(synchrony[0])
            
            light_lag[int(i), channels[ch]]=synchrony[1]*winsize
            
            
            
            
            ##stimulaion period
            #mark stim period on MA curve 
            
           
            

            ###PSTH average for each stimulus per channel
            psthal=PSTH(to, i, channels[ch], spikes_led, LED_a, LED_a)
            
            
            
            
            
            if len(succLED_a)>0:
            
            
                psthsucc=PSTH(to, i, channels[ch], spikes_led, succLED_a, LED_a)
                
                print(psthal[0].shape, psthsucc[0].shape, 'psthsshape', allpsth.shape)
            
            ###dictionary:
            
            
            
         
            
            
            psthearly=psthal[0][1:5, :] ###first 5 stimulus 
            psthlate=psthal[0][5:, :] ##rest
            binss=psthal[1] ### bins 
            winsizeled=psthal[2]  ###
            
            window=np.diff(LED_a).min()
            
            #print(window, winsizeled, 'window, winsizeled', psthlate.shape)
    
            if window>=999900:
            
                #t50ms=int(50//(1000/100))+1
                #t100ms=int(100//(1000/100))+1
                #t150ms=int(150//(1000/100))+1
                
                tpoints=gen_poin(window/1000, winsizeled)
                
                ###here window in ms and winsize is 
                
               
            else:
                
                tpoints=gen_poin(window//1000, winsizeled)
                
              
                
                
                
            psthpeaks=CatPsth(psthearly, tpoints)[0]
            
           
            
            wellpsth[channels[ch], :]=psthpeaks ##not working
            
            
            if osuc==True:
                
                if (len(succLED_a)>0):
                    
                
                    allpsth[int(i), channels[ch]]=np.mean(psthsucc[0], axis=0) 
                else: 
                    
                    allpsth[int(i), channels[ch]]=np.zeros([psthal[0].shape[1]]).tolist()
                    
            else:
                
                allpsth[int(i), channels[ch]]=np.mean(psthsal[0], axis=0) 
                
                
                    #psthpeaks
            
            binss=psthal[1]
            
            print(binss.shape, 'binssshape')
            #single channel
            

            

            ##ISI for each channel for each condition  ###LED >FB<GBJLDFIGJBDFIOJDLIGO

            isc=np.log10(np.diff(np.array(spikes_c[spikes_c<LED_s]))) 
            isc=isc[isc>=0].tolist()#pre ISI


            #extending to list one
            isic.extend(isc) #extending to list second
            gp.extend(np.repeat(i, len(isc)).tolist())
            wls.extend(np.repeat(group_iN[int(i)], len(isc)).tolist())
            light.extend(np.repeat('Off', len(isc)).tolist())

            iscled=np.log10(np.diff(np.array(spikes_c[(spikes_c>LED_s) & (spikes_c<LED_e)]))) #during whole stimulation period
            iscled=iscled[iscled>=0]

            isic.extend(iscled)
            gp.extend(np.repeat(i, len(iscled)).tolist())
            wls.extend(np.repeat(group_iN[int(i)], len(iscled)).tolist())
            light.extend(np.repeat('On', len(iscled)).tolist())

            iscledpost=np.log10(np.diff(np.array(spikes_c[spikes_c>LED_e]))) #post stimulation
            iscledpost=iscledpost[iscledpost>=0]
            
            if len(plotsi)>0:
                
                axpsth.plot(binss, np.mean(psthearly, axis=0), color=colors[ch])
                
                axpsthlast90.plot(binss, np.mean(psthlate, axis=0), color=colors[ch])
                
                #ax.set_text(((x00 + x01)/2), y+ h, "****", ha='center', va='bottom', color=col, fontsize=8)
                
                isicount, isibins=np.histogram(isc, bins=np.arange(0, 20, 1))
                isibins=[i+np.diff(isibins)[0] for i in isibins][:-1]

                
                
                logisiaxes.plot(isibins, isicount, color=colors[ch], linestyle='solid')
                
                isicount, isibins=np.histogram(iscled, bins=np.arange(0, 20, 1))
                isibins=[i+np.diff(isibins)[0] for i in isibins][:-1]
                
                logisiaxes.plot(isibins, isicount, color=colors[ch], linestyle='dashed')
                
                isicount, isibins=np.histogram(iscledpost, bins=np.arange(0, 20, 1))
                isibins=[i+np.diff(isibins)[0] for i in isibins][:-1]
                
                logisiaxes.plot(isibins, isicount, color=colors[ch], linestyle='dotted')
                
                logisiaxes.set_xlabel('Bin center', fontsize=55)
                logisiaxes.set_ylabel('Count', fontsize=55)
                logisiaxes.tick_params(axis='x', labelsize=40) 
                logisiaxes.tick_params(axis='y', labelsize=40)

                MAaxes.plot(bins_conv, color=colors[ch])

                MAaxes.set_xlabel('Bins', fontsize=55)
                MAaxes.set_ylabel('N of Spikes', fontsize=55)


                MAaxes.tick_params(axis='x', labelsize=40) 
                MAaxes.tick_params(axis='y', labelsize=40)
        
            
                

   
        arr1=arr_m_pre[int(i), :]+0.01
        arr2=arr_m_led[int(i), :]+0.01


        rch=(arr2-arr1)/arr1

        bestchid=np.argmax(rch)
        
        #axpsthlast90=WilcPSTH(wellpsth, axpsthlast90)
        
      
        
        
        presttcs=Channel_STTC(spikesynpre, start, LED_s)
        poststtcs=Channel_STTC(spikesynpost, LED_e, end)
        lightsttcs=Channel_STTC(spikesynstim, LED_s, LED_e)
        
        wellssttcpre[int(i), :]=presttcs  ##those max 12 sttc from well 
        wellssttcpost[int(i), :]=poststtcs
        wellssttclight[int(i), :]=lightsttcs
        
        
        allpsthratios=psthratios(allpsth)
        
        ###
        
       
        Rasters[indexwell, :]= rasters[bestchid, :]
        if i in LED_IDs:
            colorsR.append('blue')
            labRasters.append('light on')
        else:
            colorsR.append('black')##
            labRasters.append('light off')
        #if len(plotsi)>0:
            
            
        rasters[-1, :]=led_axis

        if i in LED_IDs:

            labrasters.append('Light on')
        else:
            labrasters.append('Light off')

        colorsrast=['black' for i in range(12)]

        #print(i, 'well')

        indexlabel=[index for index in range(len(LED_dict['Wells'])) if welllabellight in
                                                                        LED_dict['Wells'][index]]



        #print(indexlabel, welllabellight, 'colorch')


        if len(indexlabel)>0:

            colorlabel=light_label[indexlabel[0]]
        else:
            colorlabel='black'
        lightwave.extend(np.repeat(colorlabel, 12).tolist())
        colorsrast=colorsrast+[colorlabel]
        #print(colorsrast, 'color')
        #print(rasters.shape, 'rasters')
        
        
        if len(plotsi)>0:
        
            rasteraxes.eventplot(rasters, lineoffsets=lineoffset, colors=colorsrast, linelengths=linelenght, linewidths=0.3)
            rasteraxes.set_xlim(xls, xle)
            rasteraxes.set_yticklabels(labrasters)

            rasteraxes.set_xticks(np.linspace(xls, xle, 5))

            rasteraxes.set_xticklabels(np.linspace(xls//1000000, xle//1000000, 5))
            rasteraxes.set_xlabel('Time(Sec)', fontsize=55)
            rasteraxes.set_title(label+"A"+str(amps)+"Dur"+str(light_dur)+'Freq'+str(light_freq)[:5], fontsize=50)

            rasteraxes.set_ylabel('Channel', fontsize=55)
            rasteraxes.tick_params(axis='x', labelsize=40) 
            rasteraxes.tick_params(axis='y', labelsize=40)

            MAaxes.plot((np.array(kernel_mfr)/len(channels)), color=colors[-1])
            #logisiaxes.plot(isibins, np.mean(arr_isi[int(i), :, :], axis=0), color=colors[-1])
            #Databind=pd.DataFrame(acf_b[int(i), :, :])



            figrast, axrast=plt.subplots(1, 1, figsize=(16, 10))

            axrast.eventplot(rasters, lineoffsets=lineoffset, colors=colorsrast, linelengths=linelenght, linewidths=0.5)

            axrast.set_xlim(xls+19000000, xls+29000000)
            axrast.set_yticklabels(labrasters)
            axrast.set_xlabel('Time', fontsize=55)

            axrast.set_xticks(np.linspace(xls+19000000, xls+29000000, 5))

            axrast.set_xticklabels(np.linspace((xls+19000000)//1000000,  (xls+29000000)//1000000, 5))

            axrast.set_yticklabels(labrasters)

            axrast.set_title(label+"A"+str(amps)+"Dur"+str(light_dur)+'Freq'+str(light_freq)[:5], fontsize=50)



            axrast.set_ylabel('Channel', fontsize=55)
            axrast.tick_params(axis='x', labelsize=40) 
            axrast.tick_params(axis='y', labelsize=40)

                #figrast.savefig(to+'/'+'Raster'+label+compound+ID+str(i)+'.png')


    
    
     
    #colorsR=colorsR+['blue']
    #figraster, Rasteraxes=plt.subplots()
    
    #Rasteraxes.eventplot(Rasters, lineoffsets=lineoffset1, colors=colorsR, linelengths=linelenght1)
    #Rasteraxes.set_xlim(xls, xle)
    #Rasteraxes.set_yticklabels(labRasters)
    #Rasteraxes.set_xlabel('Time (microsecond)', fontsize=55)
    #Rasteraxes.set_title(label, fontsize=70)

    #Rasteraxes.set_ylabel('Channel', fontsize=55)
    #Rasteraxes.tick_params(axis='x', labelsize=40) 
    #Rasteraxes.tick_params(axis='y', labelsize=40)
    #figraster.savefig(to+'/'+'Rasters'+" "+ID+label+'.tif', format='tif', bbox_inches='tight')
    
    ##psthratios=psthratios(allptsh)
    I['ISI']=isic
    I['Well ID']=gp
    I['Group']=wls
    I['Stim']=light
    I=pd.DataFrame(I)
    
    I['Experiment']=np.repeat(ID, len(I))
    I['Compound']=np.repeat(compound, len(I))
    I['Label']=np.repeat(label, len(I))
    
    I['StimWells']=I['Well ID'].apply(lambda x : 'on' if x in LED_IDs else 'of')
    
   
    
   
    
    return (arr_m_post, arr_m_pre, I, np.array(wells), arr_m_ch, lightwave, plotsi, [LED_s, LED_e], LED_IDs, arr_m_led, (allpsth, allpsthratios, binss), (light_corr, light_lag), (wellssttcpre, wellssttcpost, wellssttclight),  arr_b[:, 1:, :]) 


# In[48]:


def gen_poin(window, winsize):
    
    
    
    #print(winsize, window, 'winsize, window')
    
    #binsize is as small as winsize value to have at least 1 bin
    #generates 50 ms windows for PSTH
    
    ### i have decided to have 5ms winsize
    
    binsize=5 #### ###change to 10ms
    
    ###and bindsize is now 50ms
    
    
    if window<50:
        
        ##if window is 50ms, binsize can be at leastn 3 tiems less. 
        
        binsize=window//3 # 6ms  
    
    
    ### if binsize is 50 and windize is 5, we have m==10
    m=int(binsize/winsize) #3 how many bins within 50ms. if 17950 is 2.6 bins 50 bins
    
    
   
    
    n50=int(window//binsize) #3 how many, 50ms bins, none, so 20 if 1000ms and 50ms binsize
    
    print(m, n50, winsize, binsize, 'm, n50, winsize, binsize')
    
    
    
    if n50<2:
        
        n50=n50+3
        
        
    
    
    
    points=[int((i*m)+1) for i in range(1, n50)]
    
    #points start in 6
    print(points, 'pointst')
    
    
    
   
    
    return points
        
        


# In[49]:


def custom_annot(pval, ax, xlocs, y, h):
    
    
    
    #print('custom plot if here')
    
    annot="{:.2f}".format(pval)
    
    x=np.mean(xlocs)
    
    #y=y+h
    
    ax.annotate('pval'+' '+annot, xy=(x, y), xycoords='data',
             xytext=(x, y+h/4), textcoords='data',
    arrowprops=dict(arrowstyle="-[",
                             linewidth=1,
                             connectionstyle="arc,armA=90,angleA=0,angleB=-40,armB=85,rad=0"),
             verticalalignment="bottom",
             horizontalalignment="left",
             fontsize=20),
             
   
        

    
    
    


# #ideas 
# 1. for the first stimulus pre 20 ms post 20 ms ratio? 
# 2. for post10 stimuli, ratio of first 20 and 20-200ms peaks if second peak exists
# 3. for post10 stimuli, ratio of first 20 and last >200ms peak if that exist
# 4. Cross-Correlation function with Light, comarison synthetic light train 
# 5. Autocorellation of stimulation window
# 6. Peaks, timepoints distribution 
# 7. Max autocorrelation distribution
# 

# In[50]:


def psthratios(psthpeaks):
    
    p20ratio=psthpeaks[:, :, 1]/psthpeaks[:, :, 0]
    p200ratio=np.zeros([p20ratio.shape[0], p20ratio.shape[1]])
    p200baseratio=np.zeros([p20ratio.shape[0], p20ratio.shape[1]])
    
    if psthpeaks.shape[2]>2:
        if np.any(psthpeaks[:, :, 2])>0:
            p200ratio=psthpeaks[:, :, 2]+0.0001/psthpeaks[:, :, 1]+0.0001
            p200baseratio=psthpeaks[:, :,  2]+0.0001/psthpeaks[:, : , 0]+0.0001

        
    return (p20ratio, p200ratio, p200baseratio)
        
    


# In[ ]:





# In[51]:


def draw_brace(ax, xspan, text):
    """Draws an annotated brace on the axes."""
    xmin, xmax = xspan
    xspan = xmax - xmin
    ax_xmin, ax_xmax = ax.get_xlim()
    xax_span = ax_xmax - ax_xmin
    ymin, ymax = ax.get_ylim()
    yspan = ymax - ymin
    resolution = int(xspan/xax_span*100)*2+1 # guaranteed uneven
    beta = 300./xax_span # the higher this is, the smaller the radius

    x = np.linspace(xmin, xmax, resolution)
    x_half = x[:resolution//2+1]
    y_half_brace = (1/(1.+np.exp(-beta*(x_half-x_half[0])))
                    + 1/(1.+np.exp(-beta*(x_half-x_half[-1]))))
    y = np.concatenate((y_half_brace, y_half_brace[-2::-1]))
    y = ymin+0.8*yspan + (.05*y - .01)*yspan # adjust vertical position

    ax.autoscale(False)
    ax.plot(x, y, color='black', lw=1)

    ax.text((xmax+xmin)/2., ymin+0.82*yspan, text, ha='center', va='bottom', fontsize=25)


# In[52]:


def draw_line(ax, xspan, text):
    """Draws an annotated brace on the axes."""
    xmin, xmax = xspan
    xspan = xmax - xmin
    ax_xmin, ax_xmax = ax.get_xlim()
    xax_span = ax_xmax - ax_xmin
    ymin, ymax = ax.get_ylim()
    yspan = ymax - ymin
    resolution = int(xspan/xax_span*100)*2+1 # guaranteed uneven
    beta = 300./xax_span # the higher this is, the smaller the radius

    x = np.linspace(xmin, xmax, resolution)
    
    
    y = ymin+0.9*yspan 
    
    y=np.repeat(y, len(x)).tolist()# adjust vertical position

    ax.autoscale(False)
    ax.plot(x, y, color='black', lw=1)

    ax.text((xmax+xmin)/2., ymin+0.9*yspan, text, ha='center', va='bottom', fontsize=18)


# In[53]:


def WilcPSTH(psthpeaks, psthax):
    
    p20=psthpeaks[:, 1] #50ms
    
    if psthpeaks.shape[1] >3:
    
        p200=psthpeaks[:, 1] #100ms
        prest=psthpeaks[:, 2] #150ms



        pvalue200=None
        pvaluerest=None
        pvalue200rest=None

        if np.any(p200)>0 and np.any(p200-p20)>0:

            res200=wilcoxon(p20, p200)



            pvalue200=res200.pvalue

            annot="{:.2f}".format(pvalue200)

            draw_brace(psthax, (20, 200), 'pval:'+annot)



        if np.any(prest)>0 and np.any(prest-p20)>0:

            resrest=wilcoxon(p20, prest)

            pvaluerest=resrest.pvalue

            #y=max(np.maximum(p20, prest))

            xmin, xmax=psthax.get_xlim()

            annot="{:.2f}".format(pvaluerest)




        if np.any(p200)>0 and np.any(prest)>0 and np.any(prest-p200)>0:

            res200r=wilcoxon(p200, prest)

            pvalue200rest=res200r.pvalue
            xmin, xmax=psthax.get_xlim()

            annot="{:.2f}".format(pvalue200rest)

            draw_brace(psthax, (200, xmax), 'pval:'+annot)



       
        
    return psthax
        
        
        
        
        
   
    
   
    
    
    
    


# In[54]:


import neo
import quantities as pq
from elephant.spike_train_correlation import spike_time_tiling_coefficient 


# In[55]:


def Channel_STTC(spikelist, start, stop):
    
    """spikelist is a list of list of channel spikes"""
    
    sttcs=[]
    
    toc=list(combinations(range(0, len(spikelist)), 2)) ##combination of indexes
    
    for tup in toc:
        
        st=STTC(spikelist[tup[0]], spikelist[tup[1]], start, stop) 
        
        sttcs.append(st)
        
    stscore=np.sqrt(np.sum([i**2 for i in sttcs]))
    
    sttcssorted=sorted(sttcs)[: : -1][:12]
    
    sttcssorted=sttcssorted+np.zeros(12-len(sttcssorted)).tolist()
    
    
    return sttcssorted
       
        
        
        
        
    
    


# In[56]:


def STTC(spike1, spike2, start, stop):
    
    ##print(stop, 'stop', spike1.dtype, spike1.shape,   spike2.shape, max(spike1), max(spike2))
    
    spiketrain1 = neo.SpikeTrain(spike1, units='us', t_stop=stop, dtype='int64')
    spiketrain2 = neo.SpikeTrain(spike2, units='us', t_stop=stop, dtype='int64')
    sttc=spike_time_tiling_coefficient(spiketrain1, spiketrain2)
    
    return sttc
    


# Channel, well, burst synchrony 
# 
# Channel sychrony: 
# STTC- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4205553/, use elephant function +
# Crosss-Correlation-Naaah
# 
# Well Synchrony-
# STTC
# 
# Network Burst-
# 
# 
# 
# 

# In[ ]:





# In[57]:


def CatPsth0(psth, tpoint):
    
    
    
    ###psth: psth with t of stim in 0 and poststimulus in 1
    ###tpoint [first peak, second peak, no peak]
    g0=psth[:, 0] #baseline
    
    g1=np.amax(psth[:, 1:tpoint[0]],  axis=1, initial=0) #first50ms
    
    l1=np.argmax(psth[:, tpoint[0]:tpoint[1]],  axis=1) ##latency of first peak
    
    
    
    g2=np.amax(psth[:, tpoint[0]:tpoint[1]], axis=1, initial=0) #f50-100ms
    
    
    g3=np.amax(psth[:, tpoint[1]:tpoint[2]], axis=1, initial=0) #100=150
    
    
    aver=[np.mean(g0, axis=0), np.mean(g1, axis=0), np.mean(g2, axis=0), np.mean(g3, axis=0)]
    #averlat=[np.mean(l1, axis=0), np.mean(l2, axis=0), np.mean(l3, axis=0)]
    
    #print(np.mean(l1, axis=0), 'latecy')
    
    return (aver, np.mean(l1, axis=0))
    
    


# In[58]:


def CatPsth(psth, tpoints):
    
    #print(tpoints, 'tpoints')
    
    
    
    ###psth: psth with t of stim in 0 and poststimulus in 1
    ###tpoint [first peak, second peak, no peak]
    
    
    g1=np.amax(psth[:, 0:tpoints[0]-1],  axis=1, initial=0)#baseline
    
    aver=[np.mean(g1, axis=0)]  ###average for the basleine? 
    l1=np.argmax(psth[:, tpoints[0]:tpoints[1]],  axis=1)
    
    for index, p in enumerate(tpoints):
        
        if index!=len(tpoints)-1:
            
            g=np.amax(psth[:, p:tpoints[index+1]], axis=1, initial=0)
            
            aver.append(np.mean(g, axis=0))
            
            
        
    return (aver, np.mean(l1, axis=0))
    
    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[1]:


def rchange(arr1, arr2):
    """Calculates the percentage change of the mean firing rate, adds 0.01 to all values to deal with zeros"""
    
    if len(arr1)==len(arr2):
    
        arr1=arr1+0.01
        arr2=arr2+0.01
    else:
        
        minN=min(len(arr1), len(arr2))
        
        arr1=arr1[:minN]
        arr2=arr2[:minN]
        
        arr1=arr1+0.01
        arr2=arr2+0.01
        



    rchange=(arr2-arr1)/arr1

    
    
    return rchange


# In[60]:


###static way of PSTH window


# In[61]:


def PSTH(to, well, channel, stim_spikes, times, timeswind):
    
    
  
    
    window=np.diff(timeswind).min()
    
    
    
    #if window>1000000:
        #window=1000000 ##look only into post-stim 1 second
   
    #window=window-10000
    
    
    winsize=5000###window//10
    
    if winsize>=5000:
        winsize=5000
        
   
    
    baselinewindow=50000 # activity before light stimulus 
    
    psths=np.zeros([len(times), int(window//winsize)+int(50000//winsize)+1])
                   ##first dimension is number of the light stimuli 
    #and 2nd dim is number of small windows
    
    if winsize<5000:
        winsize=5000
        baselinewindow=50000
        
        
        psths=np.zeros([len(times), int(window//winsize)+int(50000//winsize)+1])
        
       
    
   
   
   
    
    for time in range(0, len(times)):
        
        stim_binned=bin(stim_spikes, times[time]-(baselinewindow), times[time]+window, winsize),
        
         
        ##print(times[time]-(5*winsize), times[time]+window, winsize, '-50000', times[time], stim_binned, psths.shape)
        
        
        
        
        #pre10000
       
        
        
                
        
        
        ##print(stim_binned.shape, 'psth stim _binned shape modified')
        #stim_smoothed=np.convolve(stim_binned, np.ones(3)/3, 'same')
        #stim_smoothed=np.cumsum(stim_binned[0])/np.arange(1, len(stim_binned[0])+1, 1)
        
        stim_binned=stim_binned[0][: psths.shape[1]]
        
        
        
        if time==0:
            
            baseline=stim_binned[0]
            
        if len(stim_binned)<psths.shape[1]:
            
            stim_binned=np.concatenate((stim_binned, np.zeros(psths.shape[1]-len(stim_binned))), axis=0)
        psths[time, :]=stim_binned
        bins=np.arange(-1*baselinewindow, window, winsize)
        
        
    
    
    bins=bins/1000
    
  
    ###attepmt to make baseline very first activty before LED, only when there is very small window. 
    
    if winsize<4000:
        
        psths[:, 0]=baseline
    
    
    return (psths, bins, winsize//1000)
    


# In[62]:


def extension0(single_exp, stat, control, control2, ID, compound,group_iN, label, Param, function):
    
    
    
 
    w=stat[3].astype(int)
    
    if function=='dist-param':
        
        single_exp['Experiment'].extend(np.repeat(ID,  len(stat[0][0, w, :].ravel())))
    
        single_exp['Group'].extend(np.repeat(group_iN[w], 12))
        single_exp['Well ID'].extend(np.repeat(w, 12))
        single_exp['Stim'].extend(np.repeat('Off',  len(stat[0][0, w, :].ravel())).tolist())
        single_exp['Label'].extend(np.repeat(label,  len(stat[0][0, w, :].ravel())).tolist())
        single_exp['Compound'].extend(np.repeat(compound,  len(stat[0][0, w, :].ravel())).tolist())
        single_exp['Channel ID'].extend(stat[4][w, :].ravel())
    
        for index, p in enumerate(Param):
            
            single_exp[p].extend(stat[0][index, w, :].ravel())
            
    else:
        single_exp[Param].extend(stat[1][w, :].ravel())
        single_exp['Rchange'].extend(rchange(control, stat[1][w, :].ravel()).tolist())
        
        
   
        single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][w, :].ravel())))

        single_exp['Group'].extend(np.repeat(group_iN[w], 12))
        single_exp['Well ID'].extend(np.repeat(w, 12))
        single_exp['Stim'].extend(np.repeat('Off',  len(stat[1][w, :].ravel())).tolist())
        single_exp['Label'].extend(np.repeat(label,  len(stat[1][w, :].ravel())).tolist())
        single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][w, :].ravel())).tolist())
        single_exp['Channel ID'].extend(stat[4][w, :].ravel())
        single_exp['PSTH20'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #20
        single_exp['PSTH200'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #200 

        single_exp['PSTH20ratio'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #20tobase
        single_exp['PSTH200ratio'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #200tobase
        single_exp['STTCpre'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #20
        single_exp['STTCpost'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #200 

        single_exp['STTCStim'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
        #20tobase
        
        single_exp['Lightcorr'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
        
                                       #200 

        single_exp['LightLag'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
        single_exp['LightWave'].extend(stat[5])
        


    #single_exp['ES'].extend(np.repeat(ES,  len(stat[1][w, :].ravel())).tolist())
    
    
    if len(stat) >7:
        
        if function=='dist-param':
            single_exp['Well ID'].extend(np.repeat(w, 12))
            single_exp['Stim'].extend(np.repeat('On',  len(stat[1][0, w, :].ravel())).tolist())
            single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][0, w, :].ravel())))
            single_exp['Group'].extend(np.repeat(group_iN[w], 12))
            single_exp['Label'].extend(np.repeat(label,  len(stat[1][0, w, :].ravel())).tolist())
            single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][0, w, :].ravel())).tolist())
            single_exp['Channel ID'].extend(stat[4][w, :].ravel())

            for index, p in enumerate(Param):
                single_exp[p].extend(stat[1][index, w, :].ravel())
               
        else:
            single_exp[Param].extend(stat[9][w, :].ravel())#also stimulated
        
            single_exp['Rchange'].extend(rchange(control2, stat[1][w, :].ravel()).tolist())
       
       
            single_exp['Well ID'].extend(np.repeat(w, 12))
            single_exp['Stim'].extend(np.repeat('On',  len(stat[1][w, :].ravel())).tolist())
            single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][w, :].ravel())))
            single_exp['Group'].extend(np.repeat(group_iN[w], 12))
            single_exp['Label'].extend(np.repeat(label,  len(stat[1][w, :].ravel())).tolist())
            single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][w, :].ravel())).tolist())
            single_exp['Channel ID'].extend(stat[4][w, :].ravel())
            single_exp['PSTH20'].extend(stat[10][0][w, :, 1].ravel()) #20
            single_exp['PSTH200'].extend(stat[10][0][w, :, 2].ravel())
            
            
            
            single_exp['PSTH20ratio'].extend(stat[10][1][0][w, :].ravel()) #20tobase
            single_exp['PSTH200ratio'].extend(stat[10][1][2][w, :].ravel()) #200tobase
            
            single_exp['STTCpre'].extend(stat[12][0][w, :].ravel())#20
            single_exp['STTCpost'].extend(stat[12][1][w, :].ravel())#200 

            single_exp['STTCStim'].extend(stat[12][2][w, :].ravel())#20tobase
            single_exp['Lightcorr'].extend(stat[11][0][w, :].ravel())#200 

            single_exp['LightLag'].extend(stat[11][1][w, :].ravel())#20t
            single_exp['LightWave'].extend(stat[5])

            
            
            
            #  ###nnot sure about this, noo should work
       
    
    return single_exp
        
    


# In[63]:


def extension2(single_exp, stat, control, control2, ID, compound,group_iN, label, Param, function, ld):
    
    
    
    w=stat[3].astype(int)
    
    if function=='dist-param':
        
        single_exp['Experiment'].extend(np.repeat(ID,  len(stat[0][0, w, :].ravel())))
        single_exp['LightWave'].extend(np.repeat(0,  len(stat[0][0, w, :].ravel())).tolist())
        
    
        single_exp['Group'].extend(np.repeat(group_iN[w], 12))
        single_exp['Well ID'].extend(np.repeat(w, 12))
        single_exp['Stim'].extend(np.repeat('Off',  len(stat[0][0, w, :].ravel())).tolist())
        single_exp['Label'].extend(np.repeat(label,  len(stat[0][0, w, :].ravel())).tolist())
        single_exp['Compound'].extend(np.repeat(compound,  len(stat[0][0, w, :].ravel())).tolist())
        single_exp['Channel ID'].extend(stat[4][w, :].ravel())
    
        for index, p in enumerate(Param):
            
            single_exp[p].extend(stat[0][index, w, :].ravel())
            
    else:
        if ld==True:
            
            
            psthdict={}
            
            cols=stat[10][0].shape[2]
            
            for col in range(cols):
                
                psthdict['PSTH'+str(col*50)]=[]
            
            psthdict['Lightcorr']=[]
            psthdict['LightLag']=[]
            psthdict['LightWave']=[]
            psthdict['STTCpost']=[]
            psthdict['STTCStim']=[]
                
                
                
            if len(single_exp['Channel ID'])<1:
                
               
                
                single_exp.update(psthdict)
                
               
                
                
            for col in range(cols):
                
               
                
                
                single_exp['PSTH'+str(col*50)].extend(np.repeat(0,  len(stat[10][0][w, :, col].ravel())).tolist())
                
        single_exp[Param].extend(stat[1][w, :].ravel())
        single_exp['Rchange'].extend(rchange(control, stat[1][w, :].ravel()).tolist())
        
        
   
        single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][w, :].ravel())))

        single_exp['Group'].extend(np.repeat(group_iN[w], 12))
        single_exp['Well ID'].extend(np.repeat(w, 12))
        single_exp['Stim'].extend(np.repeat('Off',  len(stat[1][w, :].ravel())).tolist())
        single_exp['Label'].extend(np.repeat(label,  len(stat[1][w, :].ravel())).tolist())
        single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][w, :].ravel())).tolist())
        single_exp['Channel ID'].extend(stat[4][w, :].ravel())
        single_exp['LightWave'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
                                       
      

        
        single_exp['STTCpre'].extend(stat[12][0][w, :].ravel())
        
        if ld==True:
            single_exp['STTCpost'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #200 

            single_exp['STTCStim'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
            #20tobase

            single_exp['Lightcorr'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())

                                           #200 

            single_exp['LightLag'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
            single_exp['LightWave'].extend(stat[5])

        
                
        
        
            
                
        


    #single_exp['ES'].extend(np.repeat(ES,  len(stat[1][w, :].ravel())).tolist())
    
    
    if ld==True:
        
        if function=='dist-param':
            single_exp['Well ID'].extend(np.repeat(w, 12))
            single_exp['Stim'].extend(np.repeat('On',  len(stat[1][0, w, :].ravel())).tolist())
            single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][0, w, :].ravel())))
            single_exp['Group'].extend(np.repeat(group_iN[w], 12))
            single_exp['Label'].extend(np.repeat(label,  len(stat[1][0, w, :].ravel())).tolist())
            single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][0, w, :].ravel())).tolist())
            single_exp['Channel ID'].extend(stat[4][w, :].ravel())
            single_exp['LightWave'].extend(np.repeat(0,  len(stat[1][0, w, :].ravel())))

            for index, p in enumerate(Param):
                single_exp[p].extend(stat[1][index, w, :].ravel())
               
        else:
            
            psthdict={}
            
            cols=stat[10][0].shape[2]
            
            for col in range(cols):
                
                psthdict['PSTH'+str(col*50)]=[]
                
                
                
            psthdict['Lightcorr']=[]
            psthdict['LightLag']=[]
            psthdict['LightWave']=[]
            psthdict['STTCpost']=[]
            psthdict['STTCStim']=[]
                
                
                
            if len(single_exp['Channel ID'])<1:
                
              
                
                single_exp.update(psthdict)
                
                
                
            for col in range(cols):
                
                single_exp['PSTH'+str(col*50)].extend(stat[10][0][w, :, col].ravel().tolist())
                
            single_exp[Param].extend(stat[9][w, :].ravel())#also stimulated
        
            single_exp['Rchange'].extend(rchange(control2, stat[1][w, :].ravel()).tolist())
       
       
            single_exp['Well ID'].extend(np.repeat(w, 12))
            single_exp['Stim'].extend(np.repeat('On',  len(stat[1][w, :].ravel())).tolist())
            single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][w, :].ravel())))
            single_exp['Group'].extend(np.repeat(group_iN[w], 12))
            single_exp['Label'].extend(np.repeat(label,  len(stat[1][w, :].ravel())).tolist())
            single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][w, :].ravel())).tolist())
            single_exp['Channel ID'].extend(stat[4][w, :].ravel())
            
            single_exp['STTCpre'].extend(stat[12][0][w, :].ravel())#20
            single_exp['STTCpost'].extend(stat[12][1][w, :].ravel())#200 

            single_exp['STTCStim'].extend(stat[12][2][w, :].ravel())#20tobase
            single_exp['Lightcorr'].extend(stat[11][0][w, :].ravel())#200 

            single_exp['LightLag'].extend(stat[11][1][w, :].ravel())#20t
            single_exp['LightWave'].extend(stat[5])
            
            
                

            
            
            
         #  ###nnot sure about this, noo should work
   
    
    return single_exp
        
    


# In[64]:


def extension(single_exp, stat, control, control2, ID, compound,group_iN, label, Param, function, ld):
    
    
    
    w=stat[3].astype('int')
    
    print(w, 'wells')
    
    if function=='dist-param':
        
        print('nothing')
        
       
            
    else:
        if ld==True:
            
            
            psthdict={}
            
            cols=stat[10][2].tolist()
            
            
            print(stat[10][0].shape, len(cols), 'psthshape, cols')  ###psth shape hym
            
            for col in range(len(cols)-1):
                
                psthdict['PSTH'+str(cols[col])]=[]
            
            psthdict['Lightcorr']=[]
            psthdict['LightLag']=[]
            
            psthdict['STTCpost']=[]
            psthdict['STTCStim']=[]
            psthdict['MFRpost']=[]
                
                
                
            if len(single_exp['Channel ID'])<1:
                
               
                
                single_exp.update(psthdict)
                
               
                
                
            for col in range(len(cols)-1):
                
               
                
                
                single_exp['PSTH'+str(cols[col])].extend(np.repeat(0,  len(stat[10][0][w, :, col].ravel())).tolist())
                
        print(stat[4][w, :].ravel(), stat[3],  np.repeat(w, 12),  'stat4, stat3, wells, meanfiringratesldtrue')
            
                
        single_exp[Param].extend(stat[1][w, :].ravel())
        single_exp['Rchange'].extend(rchange(control, stat[1][w, :].ravel()).tolist())
        
        
   
        single_exp['Experiment'].extend(np.repeat(ID,  len(stat[1][w, :].ravel())))

        single_exp['Group'].extend(np.repeat(group_iN[w], 12))
        single_exp['Well ID'].extend(np.repeat(w, 12))
        single_exp['Stim'].extend(np.repeat('Off',  len(stat[1][w, :].ravel())).tolist())
        single_exp['Label'].extend(np.repeat(label,  len(stat[1][w, :].ravel())).tolist())
        single_exp['Compound'].extend(np.repeat(compound,  len(stat[1][w, :].ravel())).tolist())
        single_exp['Channel ID'].extend(stat[4][w, :].ravel())
        
        if ld!=True:
            
            single_exp['LightWave'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
                                       
      

        
        single_exp['STTCpre'].extend(stat[12][0][w, :].ravel())
        
        if ld==True:
            single_exp['STTCpost'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist()) #200 

            single_exp['STTCStim'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
            #20tobase

            single_exp['Lightcorr'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())

                                           #200 

            single_exp['LightLag'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())
            single_exp['LightWave'].extend(stat[5])
            single_exp['MFRpost'].extend(np.repeat(0,  len(stat[1][w, :].ravel())).tolist())

        
                
        
        
            
                
        


    #single_exp['ES'].extend(np.repeat(ES,  len(stat[1][w, :].ravel())).tolist())
    
    
    if ld==True:
        
        if function=='dist-param':
            
            print('nothingstim')
           
        else:
            
            psthdict={}
            
            cols=stat[10][2].tolist()
            print(stat[10][0].shape, len(cols), 'psthshape, cols')
            
            for col in range(len(cols)-1):
                
                psthdict['PSTH'+str(cols[col])]=[] ####here 10 is the binsize
                
                
                
            psthdict['Lightcorr']=[]
            psthdict['LightLag']=[]
            
            psthdict['STTCpost']=[]
            psthdict['STTCStim']=[]
            psthdict['MFRpost']=[]
                
                
                
            if len(single_exp['Channel ID'])<1:
                
              
                
                single_exp.update(psthdict)
                
                
                
            for col in range(len(cols)-1):
                
                single_exp['PSTH'+str(cols[col])].extend(stat[10][0][w, :, col].ravel().tolist())
                
            print(stat[4][w, :].ravel(), stat[3],  np.repeat(w, 12),  'stat4, stat3, wells, meanfiringratesldtrue')
            
            
                
            single_exp[Param].extend(stat[9][w, :].ravel())#also stimulated
        
            single_exp['Rchange'].extend(rchange(control2, stat[1][w, :].ravel()).tolist())
       
       
            single_exp['Well ID'].extend(np.repeat(w, 12))
            single_exp['Stim'].extend(np.repeat('On',  len(stat[9][w, :].ravel())).tolist())
            single_exp['Experiment'].extend(np.repeat(ID,  len(stat[9][w, :].ravel())))
            single_exp['Group'].extend(np.repeat(group_iN[w], 12))
            single_exp['Label'].extend(np.repeat(label,  len(stat[9][w, :].ravel())).tolist())
            single_exp['Compound'].extend(np.repeat(compound,  len(stat[9][w, :].ravel())).tolist())
            single_exp['Channel ID'].extend(stat[4][w, :].ravel())
            
            single_exp['STTCpre'].extend(stat[12][0][w, :].ravel())#20
            single_exp['STTCpost'].extend(stat[12][1][w, :].ravel())#200 

            single_exp['STTCStim'].extend(stat[12][2][w, :].ravel())#20tobase
            single_exp['Lightcorr'].extend(stat[11][0][w, :].ravel())#200 

            single_exp['LightLag'].extend(stat[11][1][w, :].ravel())#20t
            single_exp['LightWave'].extend(stat[5])
            single_exp['MFRpost'].extend(stat[0][w, :].ravel())
            
            
                

            
            
            
         #  ###nnot sure about this, noo should work
   
    
    return single_exp
        
    


# In[65]:


##we need to 'MFR' like calculaton 


# In[66]:


def Light_corr(binned_spikes, binned_LED):
    #autocorrelation function of the spikes
    #cross-correlation function for each channnel
    crosscorr=scipy.signal.correlate(binned_spikes, binned_LED, mode='full')
    cclags=scipy.signal.correlation_lags(binned_spikes.size, binned_LED.size, mode="full")
    ccmaxlag=cclags[np.argmax(crosscorr)]
    
    
    
    
    
    return (crosscorr, ccmaxlag)

    
    
    
    
    
    


# In[67]:


def stat_test(to, Type, df, ID, column, param, method):
    
    pairs_p=U_test(df, column, param, method, [])[1]
    pvalues=list(pairs_p.values())
    pairs=list(pairs_p.keys())
    if len(pairs)>0:
        
        formatted_pvalues = [f"p={p:.2e}" for p in pvalues]

        df=df.sort_values(by=[column])
        #vals=df['Well ID'].values
        #indexes=np.unique(vals, return_index=True)[1]
       # myorder=[vals[index] for index in sorted(indexes)]
        ploting_parameters={'x':column, 'y':param, 'data':df}


            # Create new plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_axes([0,0,1,1])

        # Plot with seaborn
        sns.boxplot(ax=ax, **ploting_parameters, showfliers=False)

        annotator = Annotator(ax, pairs, **ploting_parameters)
        #annotator.set_custom_annotations(formatted_pvalues)
        annotator.set_pvalues(pvalues)
        annotator.annotate()
        plt.savefig(to+'/'+ID+param+Type+'bygroupsstatbox'+'.tif', format='tif', bbox_inches='tight')
        plt.show()

    #plt.figure(figsize=(10, 8))
    #sns.barplot(x="Well ID", y=param, data=df)#, hue='Group', dodge=False)

    #plt.show()
    
    #sns.catplot(x="Group", y="Value", data=df, hue='Wells', kind='box', legend=False, aspect=3)

    #sns.stripplot(x="label", y="Value", data=df, color='0')
    
    #plt.savefig(to+'/'+ID+param+Type+'groupsstatbar'+'.tif', format='tif', bbox_inches='tight')

    #plt.show()
    
    return df


# In[68]:


def Burst_analysis1(to, ID, compound, series,  start, end, s, e, 
                  label_index, ch_max, label, Type,  group_iN, Param, LED_dict, ld):
    
    
    
    LED_list=[]
    LED_IDs=[]
    
    if ld==True:
        
        
        if 'Label' in LED_dict.keys():
        
            LED_IDs=[]



            # this is LED metadata 

            for lw in LED_dict['Wells']:


                LED_labels=lw

                #extracts and modifies stimulated wells, also modify this
                led_wells=label_to_id(LED_labels, WellID, Welllabel)[0]

                LED_IDs.extend(led_wells)

            light_label=LED_dict['Label']

            light_label=[string[-6:].lower() for string in light_label]

            #print(light_label, 'light_label')





            light_dur=LED_dict['Duration'][0][0]/1000


            LED=np.unique(np.array(LED_dict['Time'][0]))

            light_freq=(1000000*len(LED))/(LED[-1]-LED[0])


            if len(listamp)>0:



                amps=find_amp_weird(ID, listamp)
            else:
                amps=14

                #in the LED stim there is no amp 

        
        else:   ###electrical stimulation metadata 

            LED=LED_dict['Times'+str(label)]*50 #electrical stimulation 

            #print("LED Times", LED)# This is the stimulation already in times#change to [str(label)]
            LED_IDs=LED_dict['Wells'] #well id should be given
            LED_IDs=[int(i) for i in LED_IDs]
            amps=np.array([int(i) for i in LED_dict['Amplitudes']])   #.split(' ')])
            amps=amps[amps<0] #negative phase of biphasic stimulus
        if LED[(LED>=start) & (LED<=end)].shape[0]>0:   #the stimulation series for the given dose

            LED_s=int(LED[(LED>=start) & (LED<=end)][0])

            LED_e=int(LED[(LED>=start) & (LED<=end)][-1])

            LED_a=LED[(LED>=start) & (LED<=end)]

    else: 
        
        LED_s=end
        LED_e=end
        LED_a=[]
        
    
   
    
    'Mean of parameters, Median of Parameter, change of Parameter,classify bursts also'
    
    
    led_thing=[]
    

    Nbin=len(np.arange(0, 20, 1))-1
    
    arr_m_b=np.zeros([len(Param), 24, 12])
    arr_m_ch=np.zeros([24, 12])-1
    
    arr_led_b=np.zeros([len(Param), 24, 12])
  
    duration=(end-start)/1000000
    
    arr_w_b=np.zeros([24, 12])
    
    isic=[] 
    gp=[]
        
    wls=[]
    lightwave=[]
    
    #for each well
        
    

    
    series['Well ID']=series['Well ID'].apply(lambda x : int(x))
    wells=np.unique(series['Well ID'].sort_values(ascending=True).values)
    
    
    
    for indexwell, i in enumerate(wells):
        
        channels=np.unique(series.loc[series['Well ID']==i]['Channel Label'].sort_values(ascending=True).values)
        channelids=np.unique(series.loc[series['Well ID']==i]['Channel ID'].sort_values(ascending=True).values)
        all_spikes=len(series.loc[series['Well ID']==i])
        
        for ch in range(len(channels)):
           
            ChiD=series.loc[(series['Channel Label']==channels[ch])  & (series['Well ID']==i)]['Channel ID'].values[0]
            
            arr_m_ch[int(i), int(channels[ch])]=ChiD
            
            for index, p in enumerate(Param):
                
                
                #e.g durations of bursts in one channel, we will take mean. 
            
            
                param_c=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)][p].sort_values(ascending=True).values.astype(int)
                
                arr_m_b[index, int(i), int(ch)]=np.mean(param_c)
                
               
                
                
                
                
                if ld==True:
                    
                    param_c_pre=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)& (series['Timestamp [µs]']<LED_s)]
                    [p].sort_values(ascending=True).values.astype(int)
                    
                    
                    
                    
                
                    arr_m_b[index, int(i), int(ch)]=np.mean(param_c_pre)
                    
                    
                    param_c_led=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)& (series['Timestamp [µs]']>LED_s) & (series['Timestamp [µs]']<LED_e)]
                    [p].sort_values(ascending=True).values.astype(int)
                    
                    arr_led_b[index, int(i), int(ch)]=np.mean(param_c_led)
                    
                    
                    indexlabel=[index for index in range(len(LED_dict['Wells'])) if welllabellight in
                                                                            LED_dict['Wells'][index]]
                    if len(indexlabel)>0:
                        colorlabel=light_label[indexlabel[0]]
                    else:
                        colorlabel='black'
                else:
                    
                    colorlabel='black'
                    
                
                        
        lightwave.extend(np.repeat(colorlabel, 12).tolist())

                
                
                
                

    return  (arr_m_b, arr_led_b, led_thing, wells, arr_m_ch, lightwave, LED_list, LED_IDs)
            
    
            
    


# In[69]:


def light_burst(Welllabel, ChID, LED_dict):
    
            
    """function that returns colorofLED otherwise black, constructor of LightWave column"""
    light_label=LED_dict['Label']
        
    light_label=[string[-6:].lower() for string in light_label]

    welllabellight=Welllabel[int(ChID)]

    indexlabel=[index for index in range(len(LED_dict['Wells'])) if welllabellight in
                                                                    LED_dict['Wells'][index]]

    if len(indexlabel)>0:

        colorlabel=light_label[indexlabel[0]]
    else:
        colorlabel='black'


    return colorlabel 


# In[70]:


def Burst_analysis(to, ID, compound, dataset,  start, end, s, e, 
                  label_index, ch_max, label, Type,  group_iN, Param, LED_dict, ld, listamp):
    
    
    
    
    LED_list=[]
    LED_IDs=[]
    
    WellID=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellID.csv.npy")
    Welllabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellLabel.csv.npy")
    Channellabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\ChannelLabel.csv.npy")
    
    dataset=dataset[dataset['Dose Label']==label] 
    
    
    
    if ld==True:
        
        
        if 'Label' in LED_dict.keys():
        
            LED_IDs=[]



            # this is LED metadata 

            for lw in LED_dict['Wells']:


                LED_labels=lw

                #extracts and modifies stimulated wells, also modify this
                led_wells=label_to_id(LED_labels, WellID, Welllabel)[0]

                LED_IDs.extend(led_wells)

            light_label=LED_dict['Label']

            light_label=[string[-6:].lower() for string in light_label]

            #print(light_label, 'light_label')





            light_dur=LED_dict['Duration'][0][0]/1000


            LED=np.unique(np.array(LED_dict['Time'][0]))

            light_freq=(1000000*len(LED))/(LED[-1]-LED[0])

            dataset['LightWave']=dataset['Channel ID'].apply(lambda x: light_burst(Welllabel, x, LED_dict))







            if len(listamp)>0:



                amps=find_amp_weird(ID, listamp)
            else:
                amps=14

                #in the LED stim there is no amp 



            if LED[(LED>=start) & (LED<=end)].shape[0]>0:   #the stimulation series for the given dose

                LED_s=int(LED[(LED>=start) & (LED<=end)][0])

                LED_e=int(LED[(LED>=start) & (LED<=end)][-1])

                LED_a=LED[(LED>=start) & (LED<=end)]

            else: 

                LED_s=end
                LED_e=end
                LED_a=[]

            dataset['Stim']=dataset['Timestamp [µs]'].apply(lambda t: 'On' if (t>LED_s and t<LED_e) else 'Off')
        
    else:
        dataset['LightWave']=np.repeat('black', len(dataset))
        dataset['Stim']=np.repeat('Off', len(dataset))
    
        
        
        
    
            
   
    dataset['Group']=dataset['Well ID'].apply(lambda x: group_iN[int(x)])
    dataset['Surprise']=dataset['Surprise'].apply(lambda x: 200 if int(np.isinf(x))==1 else x)
    #dataset['Set']=dataset['File'].apply(lambda x : int(x[-5:-4]))
    dataset['Method']=dataset['File'].apply(lambda x : x[x.find('h5')+2:x.find('Bursts')])
    
    dataset['Spike Count']=np.log10(dataset['Spike Count'].values)
    
    dataset.columns = dataset.columns.str.replace("Compound ID", "Compound")
    
    dataset.columns = dataset.columns.str.replace("Dose Label", "Label")
    
    comps=np.unique(dataset['Compound'].values)
    
    
    
   
    #plt.rcParams.update({'font.size': 50})
    
    
    
    
    
   
    
    
   
    
    
    
    
    
   
    
    return (dataset, LED_list, LED_IDs)
    
    
    
    
        


# In[71]:


def Burst_analysis2(to, ID, compound, series,  start, end, s, e, 
                  label_index, ch_max, label, Type,  group_iN, Param, LED, ld):
    
    
    
    LED_list=[]
    LED_IDs=[]
    
    if ld==True:
        
        
        if 'Label' in LED_dict.keys():
        
            LED_IDs=[]



            # this is LED metadata 

            for lw in LED_dict['Wells']:


                LED_labels=lw

                #extracts and modifies stimulated wells, also modify this
                led_wells=label_to_id(LED_labels, WellID, Welllabel)[0]

                LED_IDs.extend(led_wells)

            light_label=LED_dict['Label']

            light_label=[string[-6:].lower() for string in light_label]

            #print(light_label, 'light_label')





            light_dur=LED_dict['Duration'][0][0]/1000


            LED=np.unique(np.array(LED_dict['Time'][0]))

            light_freq=(1000000*len(LED))/(LED[-1]-LED[0])


            if len(listamp)>0:



                amps=find_amp_weird(ID, listamp)
            else:
                amps=14

                #in the LED stim there is no amp 

        
        else:   ###electrical stimulation metadata 

            LED=LED_dict['Times'+str(label)]*50 #electrical stimulation 

            #print("LED Times", LED)# This is the stimulation already in times#change to [str(label)]
            LED_IDs=LED_dict['Wells'] #well id should be given
            LED_IDs=[int(i) for i in LED_IDs]
            amps=np.array([int(i) for i in LED_dict['Amplitudes']])   #.split(' ')])
            amps=amps[amps<0] #negative phase of biphasic stimulus
        if LED[(LED>=start) & (LED<=end)].shape[0]>0:   #the stimulation series for the given dose

            LED_s=int(LED[(LED>=start) & (LED<=end)][0])

            LED_e=int(LED[(LED>=start) & (LED<=end)][-1])

            LED_a=LED[(LED>=start) & (LED<=end)]

    else: 
        
        LED_s=end
        LED_e=end
        LED_a=[]
        
    
   
    
    'Mean of parameters, Median of Parameter, change of Parameter,classify bursts also'
    
    
    led_thing=[]
    

    Nbin=len(np.arange(0, 20, 1))-1
    
    ##arr_m_b=np.zeros([len(Param), 24, 12])
    
    arr_m_b={}
    arr_m_ch=np.zeros([24, 12])-1
    
    arr_led_b=np.zeros([len(Param), 24, 12])
  
    duration=(end-start)/1000000
    
    arr_w_b=np.zeros([24, 12])
    
    isic=[] 
    gp=[]
        
    wls=[]
    lightwave=[]
    
    #for each well
        
    

    
    series['Well ID']=series['Well ID'].apply(lambda x : int(x))
    wells=np.unique(series['Well ID'].sort_values(ascending=True).values)
    
    
    
    for indexwell, i in enumerate(wells):
        
        channels=np.unique(series.loc[series['Well ID']==i]['Channel Label'].sort_values(ascending=True).values)
        channelids=np.unique(series.loc[series['Well ID']==i]['Channel ID'].sort_values(ascending=True).values)
        all_spikes=len(series.loc[series['Well ID']==i])
        
        for ch in range(len(channels)):
           
            ChiD=series.loc[(series['Channel Label']==channels[ch])  & (series['Well ID']==i)]['Channel ID'].values[0]
            
            arr_m_ch[int(i), int(channels[ch])]=ChiD
            
            for index, p in enumerate(Param):
                
                
                #e.g durations of bursts in one channel, we will take mean. 
            
            
                param_c=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)][p].sort_values(ascending=True).values.astype(int)
                
                arr_m_b[index, int(i), int(ch)]=np.mean(param_c)
                
                
                
                
                if ld==True:
                    
                    param_c_pre=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)& (series['Timestamp [µs]']<LED_s)]
                    [p].sort_values(ascending=True).values.astype(int)
                    
                    
                    
                    
                
                    arr_m_b[index, int(i), int(ch)]=np.mean(param_c_pre)
                    
                    
                    param_c_led=series.loc[(series['Channel Label']==channels[ch]) & (series['Dose Label']==label)
                                    & (series['Well ID']==i)& (series['Timestamp [µs]']>LED_s) & (series['Timestamp [µs]']<LED_e)]
                    [p].sort_values(ascending=True).values.astype(int)
                    
                    arr_led_b[index, int(i), int(ch)]=np.mean(param_c_led)
                    
                    
                  
                    
                
                        
        

                
                
                
                

    return  (arr_m_b, arr_led_b, led_thing, wells, arr_m_ch, lightwave, LED_list, LED_IDs)
            
    
            
    


# In[72]:


def U_test(df, column, param, mode, pairs):
    """mode 'W'or else"""
    
    
    Exp={}
    Grp={}
    
    """If mode is W-well, each well is a singe point, if not each channel. Column is a string"""
    if mode=='W':
        df['Well ID']=df['Well ID'].apply(lambda x: str(x))
        Wells=np.unique(df['Well ID'])

        for comb in combinations(Wells, 2):

            x=df[df['Well ID']==comb[0]][param].values
            y=df[df['Well ID']==comb[1]][param].values


            if np.array(( x == y)).all() == True: 

                p=1

            else:
                stat, p=stats.ttest_rel(x, y)



            Exp[comb]=p

            if p<=0.05:

                print (comb, p)



    grp={}
    gr=np.unique(df[column].values)
    if len(gr)>0:
        
        if len(pairs)==0:
            cmbs=combinations(gr, 2)
        else:
            cmbs=pairs


        for cmb in cmbs:
            xg=df[df[column]==cmb[0]][param].values
            yg=df[df[column]==cmb[1]][param].values


            stat, p=stats.wilcoxon(xg, yg)
            
            #p=multipletests(p, alpha=0.05, method='bonferroni')[1][0]
         
            grp[cmb]=p

            if p<=0.05:

                print(cmb, p)

            
    
    
    return Exp, grp


# In[73]:


def logMFR(spikes, duration):

    mfr=(len(spikes)+1)/duration#at least one spike

    return mfr


# In[74]:


def channel_select(df, param, column, by, MC, Filter_Channel, upper):
    
    #param parameter, column for selection(channel id or well id), by (group)
    
    
    """Intragroup Channel Selection"""
    
    if Filter_Channel==True:
    #groups=np.unique(df['Group'].values)
        byg=df.groupby(by)
        pyg=df.groupby(by)#group by group, so each group filtered separately
        # lower bound ounlier removal based on 1.5IQR not  that here should be logMFR value
        if MC==False:
            
            if upper==True:
                filtered=byg.apply(lambda g:(g[(g[param]<(np.quantile(g[param].values, 0.25)-1.5*iqr(g[param].values))) & (g[param]>(np.quantile(g[param].values, 0.75)+1.5*iqr(g[param].values)))]))
            else:
                filtered=byg.apply(lambda g:(g[g[param]<(np.quantile(g[param].values, 0.25)-1.5*iqr(g[param].values))]))
            unselect=np.unique(filtered.reset_index(drop=True)[column].values)
            
        else:
            
            if upper==True:
                filtered_lower=byg.apply(lambda g:
                                   (g[g[param]<(np.quantile(g[param].values, 0.25)-(1.5*iqr(g[param].values)
                *(e**(-4*medcouple(g[param].values)))))]) if medcouple(g[param].values)>0 
                                   else (g[g[param]<(np.quantile(g[param].values, 0.25)
                -(1.5*iqr(g[param].values)*(e**(-3*medcouple(g[param].values)))))]))
                filtered_upper=pyg.apply(lambda g:
                                   (g[g[param]>(np.quantile(g[param].values, 0.75)+(1.5*iqr(g[param].values)
                *(e**(3*medcouple(g[param].values)))))]) if medcouple(g[param].values)>0 
                                   else (g[g[param]>(np.quantile(g[param].values, 0.75)
                +(1.5*iqr(g[param].values)*(e**(4*medcouple(g[param].values)))))]))
                
                unselect=np.append(np.unique(filtered_lower.reset_index(drop=True)[column].values), np.unique
                                                            (filtered_upper.reset_index(drop=True)[column].values), axis=0)
                
                #print(unselect, unselect.shape)
                
            else:
                
                filtered=byg.apply(lambda g:(g[g[param]<(np.quantile(g[param].values, 0.25)
                                                             -(1.5*iqr(g[param].values)
                    *(e**(-4*medcouple(g[param].values)))))]) 
                                 if medcouple(g[param].values)>0 else (g[g[param]<(np.quantile(g[param].values, 0.25)
                                -(1.5*iqr(g[param].values)*(e**(-3*medcouple(g[param].values)))))]))
                unselect=(filtered.reset_index(drop=True)[column].values)
        df=df[~df[column].isin(np.unique(unselect))]
    else: 
        df=df
        
    return df


# In[75]:


def weighted_SpikeNumber(df, group_iN, arsr, duration, filt):
    
    #check first channels
    
    wellid=df['Well ID'].values[0]

    channels=np.unique(df['Channel ID'].values)
    comb=combinations(channels, 2)
    rep=[]
    counts=[]
    
    
    for c in comb:
        ch1=np.unique(df[df['Channel ID']==c[0]]['Timestamp [µs]'])
        ch2=np.unique(df[df['Channel ID']==c[1]]['Timestamp [µs]'])
        ch=np.unique(df[(df['Channel ID']==c[1]) | (df['Channel ID']==c[0])]['Timestamp [µs]'])
        if len(ch)/(len(ch1)+len(ch2))<0.20:
            rep.append(sorted(list(c)))
            
    all_rep=[item for sublist in rep for item in sublist]                   
    for l in rep:
        removeflag=True
        
        for i in l:
            if i in [item for sublist in rep for item in sublist  if sublist!=l]:
                removeflag=False
        if removeFlag==True:
            rep.remove(l)
            ch=np.unique(df[(df['Channel ID']==l[1]) | (df['Channel ID']==l[0])]['Timestamp [µs]'])
            counts.append(len(ch))
           
    urep=np.unique(np.array(rep))
    allunique=np.unique(df[df['Channel ID'].isin(urep)]['Timestamp [µs]'])
    
    counts.append(len(allunique))
    
    
    for channel in channels:
        if channel not in all_rep:
            
            spks=np.unique(df[df['Channel ID']==channel]['Timestamp [µs]'])
            counts.append(len(spks))
            
    duration=duration/60 #in minute(duration is already /1000000)
    weights=[count/sum(counts) for count in counts]
    
    spike_counts=np.array(weights)*np.array(counts)
    
    if filt==True:
        
        q1=np.quantile(spike_counts, 0.25)
        q3=np.quantile(spike_counts, 0.75)
        d=1.5*iqr(spike_counts)
        upper_bound=q3+d
        lower_bound=q1-d
        
        spike_counts=[value for value in spike_counts if ((value<upper_bound) and (value>lower_bound))]
            
    arsr['awsr(min)'].append(np.sum(spike_counts)/duration) #weighted against the well
    
   
    
    
    arsr['Well ID'].append(wellid)
    arsr['Group'].append(group_iN[int(wellid)])
    
    
    return arsr


# In[76]:


def ARSR(spike_df, group_iN, duration, filt):
    
    """Array wide spiking activiity, returns Nspikes for each well"""# input is single recording initialdataframe
    arsr={'Well ID':[], 'awsr(min)':[], 'Group':[]}
    wells=spike_df.groupby('Well ID')
    
    WNspikes=wells.apply(lambda w:weighted_SpikeNumber(w.reset_index(drop=True), group_iN, arsr, duration, filt))
    return WNspikes


# In[77]:


def group_assign(df_mfawr, user_pharm, equal):
    """Experimental group assignements, input mfr/awsr dataframe for each well, and user defined nG"""
    if equal==False:
        Assignements={}
        biog=np.unique(df_mfawr['Group'].values) #Pharm assign for each group/cell line
        for b in biog:
            combin=[]

            df=df_mfawr[df_mfawr['Group']==b] 

            wells=np.unique(df['Well ID'].values)
           



            for key in user_pharm.keys(): #compounds or stimulation to be applied
                N=user_pharm[key] 
                N_auto=len(wells)//len(list(user_pharm.keys()))  #The N for the intervention
                if N>N_auto:
                    
                    N=N_auto
                   
                comb=list(combinations(wells, N))
                combin.append([list(el) for el in comb])
                #should return int
            counts=list(product(*combin))

            counts=[list(item) for item in counts if len(set(list(chain(*list(item)))))==len(list(chain(*list(item))))]
            
            fvalues=[]
            for groupz in counts:
                
                
                groupv=[df[df['Well ID'].isin(l)]['awsr(min)'].values.tolist() for l in groupz]
                

                fvalue, pvalue = stats.f_oneway(*groupv)

                fvalues.append(fvalue)
                #print(fvalues, 'fvalues')
            minindex=np.argmin(np.array(fvalues))

            assignements=counts[minindex]
            Assignements[b]=assignements
    else:
        Assignements={}
        biog=np.unique(df_mfawr['Group'].values) #Pharm assign for each group/cell line
        for b in biog:
            combin=[]

            df=df_mfawr[df_mfawr['Group']==b] 

            wells=np.unique(df['Well ID'].values)
           



            for key in user_pharm.keys(): #compounds or stimulation to be applied
               
                N=len(wells)//len(list(user_pharm.keys()))  #The N for the intervention
                
                   
                comb=list(combinations(wells, N))
                combin.append([list(el) for el in comb])
                #should return int
            counts=list(product(*combin))

            counts=[list(item) for item in counts if len(set(list(chain(*list(item)))))==len(list(chain(*list(item))))]
            
            fvalues=[]
            for groupz in counts:
                
                
                groupv=[df[df['Well ID'].isin(l)]['awsr(min)'].values.tolist() for l in groupz]
                

                fvalue, pvalue = stats.f_oneway(*groupv)

                fvalues.append(fvalue)
                #print(fvalues, 'fvalues')
            minindex=np.argmin(np.array(fvalues))

            assignements=counts[minindex]
            Assignements[b]=assignements
            
            
            
    return Assignements


# In[78]:


def match_f(df_mfr, param):
    
    df=df_mfr.sort_values(by=[param], ascending=False)
    n = len(df)-1
    matched=[]
    
    
    for i, (ind, row) in enumerate(df.iterrows()):
        #print(i, (ind, row))
    # Match the most similar untreated record to each treated record
    
        # Find the closest untreated match among records sorted 
        # higher. 'equal_or_above would' be more accurate but 
        # used 'above' for brevity        
        if i<n:
            above = df.iloc[i:]
            
            #print(above, 'above')
            wellid= df.loc[ind, 'Well ID']
            
            #print(wellid, 'wellid')
            control_above = above[(above['Well ID']!=wellid) & (~above['Channel ID'].isin(matched))]
            #print(control_above, 'control_above')
            
            if len(control_above)>0:
                match_above = control_above.iloc[0]
                #print(match_above, 'match_above')
                distance_above = match_above[param] - row[param]
                df.loc[ind, 'match'] = match_above['Well ID']
                df.loc[ind, 'distance'] = distance_above
                df.loc[ind, 'match channel'] = match_above['Channel ID']

        # Find the closest untreated match among records sorted 
        # lower. 'equal_or_below' would be more accurate but 
        # used 'below' for brevity  
        if i>0:
            below = df.iloc[:i-1]
            
            #print(below, 'below')
            wellid= df.loc[ind, 'Well ID']
            control_below = below[(below['Well ID']!=wellid)  & (~below['Channel ID'].isin(matched))]
            
            #print(control_below, 'control_below')
            
            if len(control_below)>0:
            
                match_below = control_below.iloc[-1]
                distance_below = match_above[param] - row[param]
                if i==n:
                    df.loc[ind, 'match'] = match_below['Well ID']
                    df.loc[ind, 'distance'] = distance_below
                    df.loc[ind, 'match channel'] = match_below['Channel ID']


                # Only overwrite if match_below is closer than match_above
                elif np.abs(distance_below)<np.abs(distance_above):
                    df.loc[ind, 'match'] = match_below['WellID']
                    df.loc[ind, 'distance'] = distance_below
                    df.loc[ind, 'match channel'] = match_below['Channel ID']
        matched.append(df.loc[ind, 'match channel'])
            
        matched.append(df.loc[ind, 'Channel ID'])
        
        #print(matched, 'matched')
    return df

        


# In[79]:


def match_analysis(d, param):
    
    df_mfr=d[d['Label']=='Control']
    
    
    df=df_mfr.sort_values(by=[param], ascending=False)
    n = len(df)-1
    matched=[]
    
    
    for i, (ind, row) in enumerate(df.iterrows()):
        #print(i, (ind, row))
    # Match the most similar untreated record to each treated record
    
        # Find the closest untreated match among records sorted 
        # higher. 'equal_or_above would' be more accurate but 
        # used 'above' for brevity        
        if i<n:
            above = df.iloc[i:]
            
            #print(above, 'above')
            wellid= df.loc[ind, 'Compound']
            
            #print(wellid, 'wellid')
            control_above = above[(above['Compound']!=wellid) & (~above['Channel ID'].isin(matched))]
            #print(control_above, 'control_above')
            
            if len(control_above)>0:
                match_above = control_above.iloc[0]
                #print(match_above, 'match_above')
                distance_above = match_above[param] - row[param]
                df.loc[ind, 'match'] = match_above['Well ID']
                df.loc[ind, 'distance'] = distance_above
                df.loc[ind, 'match channel'] = match_above['Channel ID']

        # Find the closest untreated match among records sorted 
        # lower. 'equal_or_below' would be more accurate but 
        # used 'below' for brevity  
        if i>0:
            below = df.iloc[:i-1]
            
            #print(below, 'below')
            wellid= df.loc[ind, 'Compound']
            control_below = below[(below['Compound']!=wellid) 
                        & (~below['Channel ID'].isin(matched))]
            
            #print(control_below, 'control_below')
            
            if len(control_below)>0:
            
                match_below = control_below.iloc[-1]
                distance_below = match_above[param] - row[param]
                if i==n:
                    df.loc[ind, 'match'] = match_below['Well ID']
                    df.loc[ind, 'distance'] = distance_below
                    df.loc[ind, 'match channel'] = match_below['Channel ID']


                # Only overwrite if match_below is closer than match_above
                elif np.abs(distance_below)<np.abs(distance_above):
                    df.loc[ind, 'match'] = match_below['WellID']
                    df.loc[ind, 'distance'] = distance_below
                    df.loc[ind, 'match channel'] = match_below['Channel ID']
        matched.append(df.loc[ind, 'match channel'])
            
        matched.append(df.loc[ind, 'Channel ID'])
    df=df.dropna()
    d=d[d['Channel ID'].isin(np.unique(df['Channel ID'].values))]
    return d
    
    
    
    


# In[80]:


def match_assign_old(matchdf):
    group1=[]
    group2=[]
    dictcount={}
    matches=matches.dropna()
    for index, i in enumerate(matches.values):
        dictcount[index]=0

        for j in matches.values:

            if i.tolist()==j.tolist() or i.tolist()==j.tolist()[::-1]:
                dictcount[index]=dictcount[index]+1


    dictcount = dict(sorted(dictcount.items(), key=lambda x: x[1], reverse=True))

    for key, value in dictcount.items():

        if matches.values[key][0] not in group2:

            group1.append(matches.values[key][0])

        if matches.values[key][1] not in group1:

            group2.append(matches.values[key][1])
    
    
    
   
    
    


# In[ ]:





# In[81]:


def match_assign(matchdf):
    grp1=[]
    grp2=[]
    for row in matches.values:



        if row[0] not in grp2:
            if row[0] not in grp1:

                if row[1] not in grp1:

                    grp1.append(row[0])
                else:
                    grp2.append(row[0])


        elif row[1]  not in grp2:
            if row[1] not in grp1:
                grp1.append(row[1])

        if row[1] not in grp1:
            if row[1] not in grp2:

                if row[0] not in grp2:

                    grp2.append(row[1])
                else:
                    grp1.append(row[1])
        elif row[0]  not in grp1:
            if row[0] not in grp2:
                grp2.append(row[0])
    return (grp1, grp2)









# In[ ]:





# In[82]:


def group_assign_II(df_mfawr, user_pharm, equal):
    """Experimental group assignements, input mfr/awsr dataframe for each well, and user defined nG, +intra f_stat"""
    
    Welllabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellLabel.csv.npy")
    WellID=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellID.csv.npy")
    #print(Welllabel)
    if equal==False:
        Assignements={}
        biog=np.unique(df_mfawr['Group'].values)
        #print(biog)#Pharm assign for each group/cell line
        for b in biog:
           

            df=df_mfawr[df_mfawr['Group']==b] 

            wells=np.unique(df['Well ID'].values)
            #print(len(wells), 'well len')
            if len(biog)>3:
                combin=[]
                for key in user_pharm.keys():#compounds or stimulation to be applied
                    N=user_pharm[key]
                    
                    #print(N, 'take')#combinations at given number 

                    comb=list(combinations(wells, N))
                    combin.append([list(el) for el in comb])
                counts=list(product(*combin))

                counts=[list(item) for item in counts if len(set(list(chain(*list(item)))))==len(list(chain(*list(item))))]
            
            else:
                combin=[]
                for i in range(100):
                    wellscurr=wells
                    comb=[]
                    
                    for key in user_pharm.keys():#compounds or stimulation to be applied
                        N=user_pharm[key]
                        #print(N, 'take')
                        'modify this as len(wells)//n of keys as if equal ns'
                        wellsind=np.random.choice(np.arange(len(wellscurr)), N, replace=False).tolist()
                        
                        comb.append(wellscurr[wellsind].tolist())
                        wellscurr=np.delete(wellscurr, wellsind)
                        
                    combin.append(comb)
                counts=combin
                #print(combin)
                
                    
            
            fvalues=[]
            for groupz in counts:
                
                
                groupv=[df[df['Well ID'].isin(l)]['awsr(min)'].values.tolist() for l in groupz]
                

                fvalue, pvalue = stats.f_oneway(*groupv)

                fvalues.append(fvalue)
        
            minindex=np.argmin(np.array(fvalues))

            assignements=counts[minindex]
            #print(assignements)
            
            
            Assignements[b]=np.unique(Welllabel[np.where(WellID==assignements)[0]])
    else:
        Assignements={}
        biog=np.unique(df_mfawr['Group'].values) #Pharm assign for each group/cell line
        for b in biog:
            combin=[]

            df=df_mfawr[df_mfawr['Group']==b] 

            wells=np.unique(df['Well ID'].values)
           



            for key in user_pharm.keys(): #compounds or stimulation to be applied
               
                N=len(wells)//len(list(user_pharm.keys()))  #The N for the intervention
                
                   
                comb=list(combinations(wells, N))
                combin.append([list(el) for el in comb])
                #should return int
            counts=list(product(*combin))

            counts=[list(item) for item in counts if len(set(list(chain(*list(item)))))==len(list(chain(*list(item))))]
            
            fvalues=[]
            for groupz in counts:
                
                
                groupv=[df[df['Well ID'].isin(l)]['awsr(min)'].values.tolist() for l in groupz]
                

                fvalue, pvalue = stats.f_oneway(*groupv)

                fvalues.append(fvalue)
               
            minindex=np.argmin(np.array(fvalues))

            assignements=counts[minindex]
            
            Assignements[b]=np.unique(Welllabel[np.where(WellID==assignements)[0]])
            
            
            
    return Assignements


# In[83]:


def group_assign_III(df_mfawr, equal):
    """Experimental group assignements, input mfr/awsr dataframe for each well, and user defined nG, +intra f_stat"""
    
    Welllabel=np.load(r"C:\Users\Amalia\Desktop\Data\Metadata\WellLabel.csv.npy")
    if equal==False:
        user_pharms=[[i, i] for i in range(4, 13)]
        fmin=[]
        nwell=[]
        for user_pharm in user_pharms:
            
            biog=np.unique(df_mfawr['Group'].values)
            #print(biog)#Pharm assign for each group/cell line
            for b in biog:


                df=df_mfawr[df_mfawr['Group']==b] 

                wells=np.unique(df['Well ID'].values)
                if len(biog)>3:
                    combin=[]
                    for key in user_pharm.keys():#compounds or stimulation to be applied
                        N=user_pharm[key] #combinations at given number 

                        comb=list(combinations(wells, N))
                        combin.append([list(el) for el in comb])
                    counts=list(product(*combin))

                    counts=[list(item) for item in counts if len(set(list(chain(*list(item)))))==len(list(chain(*list(item))))]

                else:
                    combin=[]
                    for i in range(100):
                        wellscurr=wells
                        comb=[]

                        for key in user_pharm:#compounds or stimulation to be applied
                            N=key
                            wellsind=np.random.choice(len(wellscurr), N,replace=False).tolist()

                            comb.append(wellscurr[wellsind].tolist())
                            wellscurr=np.delete(wellscurr, wellsind)

                        combin.append(comb)
                    counts=combin



            fvalues=[]

            for groupz in counts:
                groupv=[df[df['Well ID'].isin(l)]['awsr(min)'].values.tolist() for l in groupz]
                fvalue, pvalue = stats.f_oneway(*groupv)

                fvalues.append(fvalue)
               
                minindex=np.argmin(np.array(fvalues))
            
            fmin.append(min(fvalues))
            nwell.append(user_pharm[0])
        
        #print(fmin)
        #print(nwell)

        plt.figure(figsize=(10, 8))
        plt.plot(nwell, fmin)
        plt.xlabel('N of wells')
        plt.ylabel('fvalue')
        plt.show()
    else:
        Assignements={}
        biog=np.unique(df_mfawr['Group'].values) #Pharm assign for each group/cell line
        for b in biog:
            combin=[]

            df=df_mfawr[df_mfawr['Group']==b] 

            wells=np.unique(df['Well ID'].values)
           



            for key in user_pharm.keys(): #compounds or stimulation to be applied
               
                N=len(wells)//len(list(user_pharm.keys()))  #The N for the intervention
                
                   
                comb=list(combinations(wells, N))
                combin.append([list(el) for el in comb])
                #should return int
            counts=list(product(*combin))

            counts=[list(item) for item in counts if len(set(list(chain(*list(item)))))==len(list(chain(*list(item))))]
            
            fvalues=[]
            for groupz in counts:
                
                
                groupv=[df[df['Well ID'].isin(l)]['awsr(min)'].values.tolist() for l in groupz]
                

                fvalue, pvalue = stats.f_oneway(*groupv)

                fvalues.append(fvalue)
                #print(fvalues, 'fvalues')
            minindex=np.argmin(np.array(fvalues))

            assignements=counts[minindex]
            Assignements[b]=Welllabel[assignements]
            
            
            
    return fmin


# In[84]:


import numpy as np
from itertools import combinations
from itertools import product


# In[ ]:





# In[85]:


def Xhue0(df, x, hue, intra=True):
    
    """x and hue are column labels for x axis and hue"""
    
    xs=np.unique(df[x].values).astype('str')
    #print('Xs', xs)
    
  
    
    hue=np.unique(df[hue].values)
    
    if len(hue)>1:
        
        #print('hues', hue)

        pairs=list(combinations(list(product(xs, hue)), 2))
        #print('Pairs', pairs)

        pairs=[pair for pair in pairs if ((pair[0][1]==pair[1][1] and pair[0][0]!=pair[1][0]) or pair[0][0]==pair[1][0])]
        if intra==False:
            pairs=[pair for pair in pairs if pair[0][1]==pair[1][1]]

        #print('End Pairs', pairs)
    else:
        pairs=list(combinations(xs, 2))
        #print(pairs, 'end pairs')
        
    return pairs
    
    
    

    
    
    
    


# In[86]:


def Xhue(df, x, hue, intra=True, series=False):
    
    
    """x and hue are column labels for x axis and hue, 
    if age or same group different drugs, doses intra=False, for first pair and series arg True"""
    
    
    
    
    
    pairsage=[]
  
    
    hue=np.unique(df[hue].values)
    
    if series!=True:
        
        xs=np.unique(df[x].values).astype('str')
    
        if len(hue)>1:

            #print('hues', hue)

            pairs=list(combinations(list(product(xs, hue)), 2))
            #print('Pairs', pairs)

            pairs=[pair for pair in pairs if ((pair[0][1]==pair[1][1] and pair[0][0]!=pair[1][0]) or pair[0][0]==pair[1][0])]
            if intra==False:
                pairs=[pair for pair in pairs if pair[0][1]==pair[1][1]]

            #print('End Pairs', pairs)
        else:
            pairs=list(combinations(xs, 2))
            #print(pairs, 'end pairs')
    else:
        "pare groups within age and between each subsequent age"
        xs=np.unique(df[x].values).astype('str')
        
        if x=='Age':
            
            xs=np.unique(df[x].values).astype('int')
            
            xs=np.sort(xs).astype('str')
        #print('Xs', xs)
        if len(hue)>1:

            #print('hues', hue)

            pairs=list(combinations(list(product(xs, hue)), 2))
            #print('Pairs', pairs)

            pairsold=[pair for pair in pairs if (((pair[0][1]==pair[1][1] and pair[0][0]!=pair[1][0]) or 
                                              pair[0][0]==pair[1][0])) and (
                (np.where(xs==pair[0][1])[0]-np.where(xs==pair[0][0])[0])==1)]
            
            pairs=[pair for pair in pairs if (pair[0][1]!=pair[1][1]) and (pair[0][0]==pair[1][0])]
            
            #print(pairs, 'pair new, intra true')
                                             
             
            
            pairsage=list(combinations(list(product(xs, hue)), 2))
            
            if x=='Age':
                
            
                pairsage=[pair for pair in pairsage if (pair[0][1]==pair[1][1]) and (
                (np.where(xs==pair[0][0])[0]-np.where(xs==pair[1][0])[0])==-1)]
            else:
                 pairsage=[pair for pair in pairsage if (pair[0][1]==pair[1][1]) and (pair[0][0]!=pair[1][0])]
                
                
            #print(pairsage,   'pairs new, intra False')
            
        else:
            pairs=list(combinations(xs, 2))
            if x=='Age':
                
                pairs=[pair for pair in pairs if (np.where(xs==pair[0][0])[0]-np.where(xs==pair[1][0])[0])==-1]
           
                
            #print(pairs, 'end pairs')
            
                
        
        
    return (pairs, pairsage)
    
    
    

    
    
    
    


# In[87]:


def Filt_group(df, filt_dict):
    """By is column oder None, groups is groups oder all """
    for key in list(filt_dict.keys()):
        
        #print('key', key)
        if key!='None':
            
            #print('key', key, 'doing')
            df=df[~df[key].isin(filt_dict[key])]
        else:
            df=df
            
    return df
    


# In[88]:


def scatter_doses(dfs, experiment, Param, to, order, colordict, plotparam):
    '''df if experiment compound selected dataframe'''
    
    from statsmodels.formula.api import ols
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    compounds=np.unique(dfs['Compound'].values)
    
    comps=[comp for comp in compounds if comp!='No Group']
    
    #print(comps)
    
    max_value=dfs[Param].mean()+5*(np.std(dfs[Param].values))
    

    
    colormap=['blue', 'grey', 'green', 'yellow', 'red']
    
    plt.style.use('ggplot')
    
    
    
   
    Labels=np.unique(dfs['Label'].values)
    #markerlist=list(Line2D.markers.keys())[:24]
    groups=np.unique(dfs['Group'].values)
    
    
    
    
    combin=list(combinations(Labels, 2))
    
    axs=plt.subplots(len(groups), len(combin), figsize=plotparam['figsize'],  sharey=plotparam['sharey'])
    
    for cindex, comb in enumerate(combin):
        
        print(cindex)
        
       
        
        
        if np.where(order==comb[0])[0]>np.where(order==comb[1])[0]:
            comb=comb[::-1]
        

        

        
        #dfdoses=dfs[(dfs['Label']==comb[0]) | (dfs['Label']==comb[1])]
        
        #sns.scatterplot(data=dfdoses, x=comb[0], y=comb[1], hue='Compound')
        
        for index, comp in enumerate(comps):
            
            
            for indg, g in enumerate(groups):


                df=dfs[(dfs['Compound']==comp) & (dfs['Group']==g)]



                x=df[df['Label']==comb[0]].sort_values(by=['Channel ID'])

                colorlist=np.repeat(colordict[g], len(x))




                #markerswell=x['Well ID'].values.tolist()
                #markers=[markerlist[int(i)] for i in markerswell]

                #colormap=np.arange(0, len(x['Channel ID'].values), 1)



                x=x[Param].values
                y=df[df['Label']==comb[1]].sort_values(by=['Channel ID'])
                y=y[Param].values



                axs[indg, cindex].scatter(x, y, color=colormap[index], label=g)
                axs[indg, cindex].set_xlabel(comb[0], fontsize=plotparam['labelsize'])
                axs[indg, cindex].set_ylabel(comb[1],  fontsize=plotparam['labelsize'])
                #plt.xlim(0, max_value)
                #plt.ylim(0, max_value)





                model = LinearRegression().fit(x.reshape(-1, 1), y)#sm.OLS(y, sm.add_constant(x)).fit()
                y_hat=model.predict(x.reshape(-1, 1))

                r2=r2_score(y, y_hat)   

                #plt.plot(x, y_hat, color=colormap[index], label="r\u00b2 = {:.3f}".format(r2)+" "+comp)
                axs[indg, cindex].plot(x, x, color='black', label='Unity', linestyle='dotted')


                #plt.annotate("r-squared = {:.3f}".format(r2), (0, 1), color=colormap[index])


                #plt.legend(fontsize=15)
                #plt.savefig(to+'/'+experiment+'Scatter'+comb[0]+comb[1]+comp+g+'.png')

                #influence = model.get_influence()
                #inf_sum = influence.summary_frame()


                #student_resid = influence.resid_studentized_external

       
        

        
        
       
        plt.savefig(to+'/'+compound+experiment+'Scatter'+'.png')
        
    return Labels


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[89]:


def stat_plot_annot(pairs, gdf, x, param, hue, to, tosub, plottype, plotTitle, hue_order_all, ageplot):
    
    """Applying statistical test, U-shape for independent and Wilcxon for dependent sampels are routinely applied
    pair argument controls this, pairs[0] are indepedent samples and pairs[1] are depenedent samples"""
    
    #print(gdf[:3], 'awsrerroe')
    
    ide=gdf['Experiment'].values[0]
    lightw=gdf['LightWave'].values[0]
    
    #print(ide, 'Experiment')
    
    hue_order=hue_order_all[ide]
    
    if plotTitle=='Group':
        Title=str(gdf['Group'].values[0])
    elif plotTitle=='Experiment':
        Title=str(gdf['Experiment'].values[0])
    else:
        Title=str(gdf['Label'].values[0])
        
    comp=np.unique(gdf['Compound'].values).tolist()
    #print(comp, 'compounds')
    
    Scomp=" ".join(comp)
    test2='Mann-Whitney'
    
    #print('X', x)
    #print('param', param)
    #print('hue', hue)
    #print('pairs[0]', pairs[0])
    #print('df', gdf[:10])
    #print('pairs[1]', pairs[1])
    
    if gdf[x].values[0].isdigit()==True:
        order=np.sort(np.unique(gdf[x].values).astype('int')).astype('str')
    elif x=='Compound':
        if 'Vehicle' in np.unique(gdf[x].values):
            ordlist=np.unique(gdf[x].values).tolist()
            ordlist.remove('Vehicle')
            order=['Vehicle'].extend(ordlist)
        else:
            order=None
    else:
        order=None
        
        

    
    if (x=='Age') and (hue=='Group') and (len(pairs[1])>1):
        free=pairs[0]
        depends=pairs[1]
        
    elif (x=='Compound') and (hue=='Label') and (len(pairs[1])>1):
        free=pairs[1]
        
       
        depends=pairs[0]
        
    ##single compound many series    
    elif (x=='Compound') and (hue=='Label') and (len(pairs[1])<1) and (len(pairs[0])>1):
        
        free=[]
        
        depends=pairs[0]
        
        
    else:
        free=pairs[0]
        ageplot=False
        depends=[]
   
        
        
    
   
    
    if (len(pairs[0])>0) or (len(pairs[1]))>0:
        
        
        if plottype=='box':
            
            fig=plt.figure(figsize=(10, 8))
            
            if hue=='Label':
                
                if len(comp)<2:
                    test2='Wilcoxon'
                
                    
                
                hue_order=hue_order[comp[0]]
                hue_order=[i[0] for i in hue_order]
                
                if param=='Change':
                    
                    hue_order=np.delete(np.array(hue_order), int(np.where(np.array(hue_order)=='Control')[0]))
                    
                
                
                #if 'Vehicle' in np.unique(gdf['Label'].values):
                    #doses=np.unique(gdf['Label'].values).remove('Vehicle', '0').tolist()
                    #doses=sorted([int(''.join(c for c in x if c.isdigit())) for x in doses])
                    #hueorder=['0', 'Vehicle'].extend([str(i) for i in doses])
                #else:
                    #doses=np.unique(gdf['Label'].values).tolist()
                    #doses=sorted([int(''.join(c for c in x if c.isdigit())) for x in doses])
                    #hueorder=['0', 'Vehicle'].extend([str(i) for i in doses])
                                                                
    
                
                y_values=gdf[param].values
            
            
            
                #print(np.unique(gdf[hue].values), 'hues', '????')
                if len(np.unique(gdf[hue].values))>1:
                    
                    
                    ax = sns.boxplot(data=gdf, x=x, y=param, hue=hue, hue_order=hue_order, order=order, showfliers = False)
                                     
                   #ax.set_ylim(np.min(y_values), np.max(y_values))
                    adjust_box_widths(fig, 0.9)
                    
                    if len(free)>0:
                        
                    
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, hue=hue, hue_order=hue_order, order=order)
                    
                        
                    
                    if ageplot==True:
                        
                        try:
                            #print('was here')
                            
                            annot2=Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue, hue_order=hue_order, 
                                              order=order)
                            annot2.configure(test='Wilcoxon', verbose=2, loc="inside", format_text='simple')

                            annot2.apply_test()
                            annot2.annotate()
                        except:
                            
                            try:
                                annot2 =Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue, 
                                                  hue_order=hue_order, order=order)
                                annot2.configure(test='Kruskal', verbose=2, loc="inside")
                                annot2.apply_test()
                                annot2.annotate()
                            except:
                                 annot2 =Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue, 
                                                  hue_order=hue_order, order=order)
                                

                        
                        
                        
                    
                    
                else:
                    ax = sns.boxplot(data=gdf, x=x, y=param, showfliers = False)
                   #ax.set_ylim(np.min(y_values), np.max(y_values))
                    adjust_box_widths(fig, 0.9)
                    
                    if len(free)>0:
                        
                        annot = Annotator(ax, free, data=gdf, x=x, y=param)
                   
                    
            else:
                
                y_values=gdf[param].values
                
                if len(np.unique(gdf[hue].values))>1:
                    ax = sns.boxplot(data=gdf, x=x, y=param, hue=hue, order=order, showfliers = False)
                   #ax.set_ylim(np.min(y_values), np.max(y_values))
                    adjust_box_widths(fig, 0.9)
                    
                    if len(free)>0:
                    
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, hue=hue, order=order)
                    
                    if ageplot==True:
                        
                    
                        try:
                            
                            annot2 = Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue, order=order)
                            annot2.configure(test='Wilcoxon', verbose=2, loc="inside", paired=False)
                            annot2.apply_test()
                            annot2.annotate()
                            
                        except:
                            
                            try:
                            
                       
                                annot2 =Annotator(ax, depends, data=gdf, x=x, y=param, 
                                                  hue=hue,hue_order=hue_order, order=order)
                                annot2.configure(test='Kruskal', verbose=2, loc="inside", text_format='simple')
                                annot2.apply_test(pairs=False)
                                annot2.annotate()
                            except:
                                annot2 =Annotator(ax, depends, data=gdf, x=x, y=param, 
                                                  hue=hue,hue_order=hue_order, order=order,showfliers=False)
                                
                                
                            

                        
                else:
                    ax = sns.boxplot(data=gdf, x=x, y=param, order=order, showfliers = False)
                    #ax.set_ylim(np.min(y_values), np.max(y_values))
                    adjust_box_widths(fig, 0.9)
                    
                    if len(free)>0:
                        
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, order=order, showfliers=False )
                    
                
            if len(pairs[0])<10:
                
                if len(free)>0:
                    try:
                        
                        

                        annot.configure(test=test2, verbose=2, loc="inside",
                                        text_format='simple')
                        annot.apply_test()
                        annot.annotate()
                    except:
                        print('No test possible')

               

                    
                
            plt.title("by"+Title+x+hue+param)
            plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
            plt.savefig(tosub+'/'+param+hue+x+Title+ide+str(lightw)+Scomp+'.png', dpi=300, bbox_inches='tight')
            
            if len(free)>0:
                annot.reset_configuration()
            
            

        if plottype=='bar':
            
            plt.figure(figsize=(10, 8))
            if hue=='Label':
                
                hue_order=hue_order[comp[0]]
                hue_order=[i[0] for i in hue_order]
    
                
                y_values=gdf[param].values
                
                if len(np.unique(gdf[hue].values))>1:
                    ax = sns.barplot(data=gdf, x=x, y=param, hue=hue, order=order, hue_order=hue_order)
                    #ax.set_ylim(np.min(y_values), np.max(y_values))
                    
                    if len(free)>0:
                    
                    
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, hue=hue, order=order, hue_order=hue_order)
                    
                    if ageplot==True:
                        
                        try:
                            annot2 = Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue, order=order, hue_order=hue_order)
                            annot2.configure(test='Wilcoxon', verbose=2, loc="inside", paired=False)
                            annot2.apply_test()
                            annot2.annotate()
                                
                        except:
                            
                            annot2 =Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue,hue_order=hue_order, order=order)
                            annot2.configure(test='Kruskal', verbose=2, loc="inside")
                            annot2.apply_test()
                            annot2.annotate()


                else:
                    ax = sns.barplot(data=gdf, x=x, y=param, order=order)
                    #ax.set_ylim(np.min(y_values), np.max(y_values))
                    
                    if len(free)>0:
                    
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, order=order)
                
                
            else:
                
                
                if len(np.unique(gdf[hue].values))>1:
                    
                    #print('order', order)
                    print('hue', hue)
                    
                    print(free, 'free')
                    
                    print(depends, 'depends')
                    
                    print(hue_order, 'hue_order')
                    
                    
                    ax = sns.barplot(data=gdf, x=x, y=param, hue=hue, order=order)
                    y_values=gdf[param].values
                    #ax.set_ylim(np.min(y_values), np.max(y_values))
                    
                    if len(free)>0:
                        
                    
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, hue=hue, order=order)
                    
                    if ageplot==True:
                        
                        
                        try:
                            
                            
                            annot2 = Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue, order=order)
                            annot2.configure(test='Wilcoxon', verbose=2, loc="inside", text_format='simple')
                            annot2.apply_test()
                            annot2.annotate()

                        except:
                            
                            try:

                            
                                annot2=Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue,
                                                 order=order)
                                annot2.configure(test='Kruskal', verbose=2, loc="inside", text_format='simple')
                                annot2.apply_test()
                                annot2.annotate()
                            except:
                                 annot2=Annotator(ax, depends, data=gdf, x=x, y=param, hue=hue,
                                                  order=order)
                                


                    
                    
                    
                else:
                    ax = sns.barplot(data=gdf, x=x, y=param, order=order)
                    y_values=gdf[param].values
                    #ax.set_ylim(np.min(y_values), np.max(y_values))
                    
                    if len(free)>0:
                        annot = Annotator(ax, free, data=gdf, x=x, y=param, order=order)
                
                
                
            if len(pairs[0])<16:
                try:
                    annot.configure(test='Kruskal', verbose=2, loc='inside', text_format='simple')
                    annot.apply_test()
                    annot.annotate()
                except:
                    'No test'

           
                
          
            
            plt.title("by"+Title+x+hue+param)
            plt.xlabel(x)
            plt.ylabel(param)
            plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))
            plt.savefig(tosub+'/'+param+hue+x+Title+ide+str(lightw)+Scomp+'.png', dpi=300, bbox_inches='tight')
            
            if len(free)>0:
                
                annot.reset_configuration()
            
    
        
    
    return pairs[0]


# In[90]:


def Plot0(dfs, to, filt_dict, effects, hue_order_all):
    folderISI=to+'/'+'ISI'
    folderMFR=to+'/'+'MFR'
    folderawsr=to+'/'+'AWSR'
    
    for directory in [folderISI, folderMFR, folderawsr]:
        if not os.path.exists(directory):
            os.makedirs(directory)

    
    
    if 'Age' in effects:
        

        for df in dfs:
            
            df=df[df['Compound']=='No Group']
            
            if len(df)>0:
            
                df=Filt_group(df, filt_dict)

                if 'ISI' in df.columns:

                    tosub=folderISI  #tobe defined 

                    pairs=Xhue(df, 'Age', 'Group', intra=False)
                    pairs=stat_plot_annot(pairs, df, 'Age', 'ISI', 'Group', to, tosub, 'box', 'Group', hue_order_all)

                if 'MFR(Hz)' in df.columns:

                    tosub=folderMFR

                    pairs=Xhue(df, 'Age', 'Group', intra=False)

                    pairs=stat_plot_annot(pairs, df, 'Age', 'MFR(Hz)', 'Group', to, tosub, 'bar', 'Group', hue_order_all)

                if 'awsr(min)' in df.columns:

                    tosub=folderawsr

                    pairs=Xhue(df, 'Age', 'Group', intra=False)

                    pairs=stat_plot_annot(pairs, df, 'Age', 'awsr(min)', 'Group', to, tosub, 'box', 'Group', hue_order_all)
                    
                    df['Age']=df['Age'].values.astype('int')
                
                    
                    if len(np.unique(df['Group'].values))>1:

                        sns.catplot(data=df, x='Age', y='awsr(min)', hue='Group', kind='point')
                    else:
                        
                        #print(df[:10])
                        #print(df['Age'].dtype)
                        sns.catplot(data=df, x='Age', y='awsr(min)', kind='point')

                    plt.savefig(tosub+'/'+'awsr(min)''Group'+'Age'+'.png', dpi=300, bbox_inches='tight')

    if 'Drug'  in effects:

        for df in dfs:
            
            #print(df.columns, 'columnsfindmfr')
            
            
            
            dfdrug=df#Filt_group(df, filt_dict)
            

            compounds=np.unique(df['Compound'].values)

            for compound in compounds:
                
                dfcompound=dfdrug[(dfdrug['Compound']==compound)]
                
                #print(compound, 'compoundformfr')
                
                if compound!='Vehicle':
                    
                   

                   
                    dfdrug=dfdrug[(dfdrug['Compound']==compound)| (df['Compound'].str.contains('Vehicle'))]
                    
                    

                    resdf=dfdrug.groupby('Experiment')
                    resdfcomp=dfcompound.groupby('Experiment')
                    
                                          

                    if 'ISI' in df.columns:
                        tosub=folderISI

                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Group', intra=True), g, 'Compound', 'ISI', 'Group', to, tosub, 'box', 'Label', hue_order_all)))

                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True), g, 'Compound', 'ISI', 'Label', to, tosub, 'box', 'Group', hue_order_all)))
                    if 'MFR(Hz)' in df.columns:

                        tosub=folderMFR
                        #print('MFRplotting')

                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True), g, 'Compound', 'MFR(Hz)', 'Label', to, tosub, 'box', 'Group', hue_order_all)))
                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Group', intra=True), g, 'Compound', 'MFR(Hz)', 'Group', to, tosub, 'box', 'Label', hue_order_all)))
                        res=resdfcomp.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True), g, 'Compound', 'MFR(Hz)', 'Label', to, tosub, 'box', 'Group', hue_order_all)))
                    if 'awsr(min)' in df.columns:

                        tosub=folderawsr
                        #print(df[:3], 'awsr')
                        res=resdfcomp.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True), g, 'Compound', 'awsr(min)', 'Label', to, tosub, 'bar', 'Group', hue_order_all)))
                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=False), g, 'Compound', 'awsr(min)', 'Label', to, tosub, 'bar', 'Group', hue_order_all)))
                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Group', intra=True), g, 'Compound', 'awsr(min)', 'Group', to, tosub, 'bar', 'Label', hue_order_all)))



    if "Cell Line" in effects:

        for df in dfs:
            
            df=Filt_group(df, filt_dict)

            df=df[df['Compound']=='No Group']
            
            if len(df)>0:

                resdf=df.groupby('Experiment')

                if 'ISI' in df.columns:
                    tosub=folderISI
                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True), g, 'Group', 'ISI', 'Well ID', to, tosub, 'box', 'Group', hue_order_all))

                if 'MFR(Hz)' in df.columns:

                    tosub=folderMFR

                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True), g, 'Group', 'MFR(Hz)', 'Well ID', to, tosub, 'box', 'Group', hue_order_all))


                if 'awsr(min)' in df.columns:
                    
                    #print(df[:3], 'awsr')

                    tosub=folderawsr
                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True), g, 'Group', 'awsr(min)', 'Well ID', to, tosub, 'box', 'Group', hue_order_all))




        
    return 'Done'
         

                 
    


# In[91]:


def Plot1(dfs, to, filt_dict, effects, hue_order_all, ageplot, Parameter):
    folderISI=to+'/'+'ISI'
    folderMFR=to+'/'+Parameter
    folderawsr=to+'/'+'AWSR'
    
    for directory in [folderISI, folderMFR, folderawsr]:
        if not os.path.exists(directory):
            os.makedirs(directory)

    
    
    if 'Age' in effects:   
        

        for df in dfs:
            
            df=df[df['Compound']=='No Group']
            
            #print(np.unique(df['Well ID'].values), 'wellfiltered?age')
            
            if len(df)>0:
            
                df=Filt_group(df, filt_dict)

                if 'ISI' in df.columns:

                    tosub=folderISI  #tobe defined 

                    pairs=Xhue(df, 'Age', 'Group', intra=False)
                    pairs=stat_plot_annot(pairs, df, 'Age', 'ISI', 'Group', to, tosub, 'box', 'Group', hue_order_all, ageplot)

                if Parameter in df.columns:

                    tosub=folderMFR

                    pairs=Xhue(df, 'Age', 'Group', intra=True, series=True)

                    pairs=stat_plot_annot(pairs,
                            df, 'Age', Parameter, 'Group', to, tosub, 'bar', 'Group', hue_order_all, ageplot)

                if 'awsr(min)' in df.columns:

                    tosub=folderawsr

                    pairs=Xhue(df, 'Age', 'Group', intra=False)

                    pairs=stat_plot_annot(pairs, df, 'Age', 'awsr(min)', 'Group', to, tosub, 'box', 'Group', hue_order_all, ageplot)
                    
                    df['Age']=df['Age'].values.astype('int')
                
                    
                    if len(np.unique(df['Group'].values))>1:

                        sns.catplot(data=df, x='Age', y='awsr(min)', hue='Group', kind='point')
                    else:
                        
                        #print(df[:10])
                        #print(df['Age'].dtype)
                        sns.catplot(data=df, x='Age', y='awsr(min)', kind='point')

                    plt.savefig(tosub+'/'+'awsr(min)''Group'+'Age'+'.png', dpi=300, bbox_inches='tight')

    if 'Drug'  in effects:

        for df in dfs:
            
            
            
            dfdrug=Filt_group(df, filt_dict)
            
           

            compounds=np.unique(df['Compound'].values)

            for compound in compounds:
                
                dfcompound=dfdrug[(dfdrug['Compound']==compound)]
                
                #print(compound, 'compoundformfr')
                
                if compound!='Vehicle':
                    
                   

                   
                    dfdrug=dfdrug[(dfdrug['Compound']==compound)| (dfdrug['Compound'].str.contains('Vehicle'))]
                    
                    

                    resdf=dfdrug.groupby(['Experiment'])
                    resdfcomp=dfcompound.groupby('Experiment')
                    
                                          

                    if 'ISI' in df.columns:
                        tosub=folderISI

                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue
                                                                            (g, 'Compound', 'Group', intra=True, series=False), g, 
                                                        'Compound', 'ISI', 'Group', to, tosub, 'box', 'Label', hue_order_all, ageplot)))

                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue
                                                                                (g, 'Compound', 'Label', intra=True), g, 'Compound', 'ISI', 'Label', to, tosub, 'box', 'Group', hue_order_all, ageplot)))
                    if Parameter in df.columns:

                        tosub=folderMFR
                        #print('MFRplotting')
                        
                        
                        
                        fg=resdf.apply(lambda g : relative(g, tosub, hue_order_all, Parameter, 'Group', compound, to, 'True'))
                        
                        
                        gf=resdf.apply(lambda g : relative(g, tosub, hue_order_all, Parameter, 'Label', compound, to, 'True'))
                        
                  
                        
                        
                        
                       
                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot
                                                                           (Xhue(g, 'Compound', 'Label', intra=True, 
                                                                                 series=True), g, 'Compound', 
                                                Parameter, 'Label', to, tosub, 'box', 'Group', hue_order_all, ageplot)))
                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot
                                                                           (Xhue(g, 'Compound', 'Group',
                                                                                 intra=True, series=False), g, 'Compound', Parameter, 'Group', to, tosub, 'box', 'Label', hue_order_all, ageplot)))
                        #print('printing?')
                        res=resdfcomp.apply(lambda x: x.groupby('Group').apply(lambda g: 
                                    stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True, series=True),
                                    g, 'Compound', Parameter, 'Label', to, tosub, 'box', 'Group', hue_order_all, True)))
                        
                        
                        
                    if 'awsr(min)' in df.columns:

                        tosub=folderawsr
                        #print(df[:3], 'awsr')
                        res=resdfcomp.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True), g, 'Compound', 'awsr(min)', 'Label', to, tosub, 'bar', 'Group', hue_order_all, ageplot)))
                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot
                                                                           (Xhue(g, 'Compound', 'Label', intra=True), g, 
                                                                            'Compound', 'awsr(min)', 'Label', to, tosub, 'bar', 'Group', hue_order_all, ageplot)))
                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound',
                                                                                                          'Group', intra=True, series=False), g, 'Compound', 'awsr(min)', 'Group', to, tosub, 'bar', 'Label', hue_order_all, ageplot)))



    if "Cell Line" in effects:
        
        

        for df in dfs:
            
            #df=Filt_group(df, filt_dict)
            
            dfpharm=df

            df=df[df['Compound']=='No Group']
            
            if len(df)>0:

                resdf=df.groupby('Experiment')

                if 'ISI' in df.columns:
                    tosub=folderISI
                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True, series=False), g, 'Group', 'ISI', 'Well ID', to, tosub, 'box', 'Group', hue_order_all, ageplot))

                if Parameter in df.columns:
                    
                    #print(np.unique(df['Well ID'].values), 'wellfiltered?')
                    
                    

                    tosub=folderMFR
                    
                    
                    
                    pairs=Xhue(df, 'Age', 'Group', intra=True, series=False)
                    
                    #print(pairs, 'pairs, what')

                    pairs=stat_plot_annot(pairs,
                            df, 'Age', Parameter, 'Group', to, tosub, 'box', 'Group', hue_order_all, False)

                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True, series=False), 
                                                              
                                                              g, 'Group', Parameter, 'Well ID', to, tosub, 'box', 
                                                              'Group', hue_order_all, ageplot))


                if 'awsr(min)' in df.columns:
                    
                    #print(df[:3], 'awsr')

                    tosub=folderawsr
                    #res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True),
                                                              #g, 'Group', 'awsr(min)', 'Well ID', to, tosub, 'box', 'Group', hue_order_all, ageplot))




        
    return 'Done'


# In[92]:


def statistical_test(df, column, param, mode, test, ns,  pairs):
    """mode 'W'or else"""
    
    
    Exp={}
    Grp={}
    
    sign=[]
    
    """If mode is W-well, each well is a singe point, if not each channel. Column is a string"""
    if mode=='W':
        df['Well ID']=df['Well ID'].apply(lambda x: str(x))
        Wells=np.unique(df['Well ID'])

        for comb in combinations(Wells, 2):

            x=df[df['Well ID']==comb[0]][param].values
            y=df[df['Well ID']==comb[1]][param].values


            if np.array(( x == y)).all() == True: 

                p=1

            else:
                stat, p=stats.ttest_rel(x, y)



            Exp[comb]=p

            #if p<=0.05:

                ##print (comb, p)



    grp={}
    gr=np.unique(df[column].values)
    if len(gr)>0:
        
        if len(pairs)==0:
            cmbs=combinations(gr, 2)
        else:
            cmbs=pairs


        for cmb in cmbs:
            xg=df[df[column]==cmb[0]][param].values
            yg=df[df[column]==cmb[1]][param].values


            stat, p=test(xg, yg)
            
            p=multipletests(p, alpha=0.05, method='bonferroni')[1][0]
         
            if ns==False:

                if p<=0.05:

                    #print(cmb, p)
                
                    grp[cmb]=p
            else:
                grp[cmb]=p
                
            if p<=0.05:
                
                sign.append(cmb)
                
                

            
    
    
    return Exp, grp, sign


# In[93]:


def Plot(dfs, to, filt_dict, effects, hue_order_all, ageplot, Parameter, Li):
    folderISI=to+'/'+'ISI'
    folderMFR=to+'/'+Parameter
    folderawsr=to+'/'+'AWSR'
    
    for directory in [folderISI, folderMFR, folderawsr]:
        if not os.path.exists(directory):
            os.makedirs(directory)

    
    
    if 'Age' in effects:   
        

        for df in dfs:
            
            #df=df[df['Compound']=='No Group']
            
            df=df[df['Label']=='Control']
            
            
            
            #print(np.unique(df['Well ID'].values), 'wellfiltered?age')
            
            if len(df)>0:
            
                df=Filt_group(df, filt_dict)

                if 'ISI' in df.columns:

                    tosub=folderISI  #tobe defined 

                    pairs=Xhue(df, 'Age', 'Group', intra=False)
                    pairs=stat_plot_annot(pairs, df, 'Age', 'ISI', 'Group', to, tosub, 'box', 'Group', hue_order_all, ageplot)

                if Parameter in df.columns:

                    tosub=folderMFR

                    pairs=Xhue(df, 'Age', 'Group', intra=True, series=True)
                    
                    #print(df.columns, 'troubleshootingcolumns')

                    pairs=df.groupby('LightWave').apply(lambda d: stat_plot_annot(pairs,
                            d, 'Age', Parameter, 'Group', to, tosub, 'box', 'Group', hue_order_all, ageplot))
                    
                                   

                if 'awsr(min)' in df.columns:

                    tosub=folderawsr

                    pairs=Xhue(df, 'Age', 'Group', intra=False)

                    pairs=stat_plot_annot(pairs, df, 'Age', 'awsr(min)', 'Group', to, tosub, 'box', 'Group', hue_order_all, ageplot)
                    
                    df['Age']=df['Age'].values.astype('int')
                
                    
                    if len(np.unique(df['Group'].values))>1:

                        sns.catplot(data=df, x='Age', y='awsr(min)', hue='Group', kind='point')
                    else:
                        
                        #print(df[:10])
                        #print(df['Age'].dtype)
                        sns.catplot(data=df, x='Age', y='awsr(min)', kind='point')

                    plt.savefig(tosub+'/'+'awsr(min)''Group'+'Age'+'.png', dpi=300, bbox_inches='tight')

    if 'Drug'  in effects:

        for df in dfs:
            
            
            
            
            
            dfdrug=Filt_group(df, filt_dict)
            
            
            dfdrugVehi=dfdrug.copy()
            
           
            
           
                                
            dfdrugVehi['Compound']=dfdrugVehi['Well ID'].apply(lambda x: dfdrugVehi[
                                     (dfdrugVehi['Well ID']==x)
                                                 & (dfdrugVehi['Compound']!='No Group')]['Compound'].values[0])
           

            compounds=np.unique(df['Compound'].values)

            for compound in compounds:
                
                dfcompound=dfdrug[(dfdrug['Compound']==compound)]
                
                #print(compound, 'compoundformfr')
                
                if ('Vehicle' not in compound) and (compound!='No Group'):
                    
                   

                   
                    dfdrugVeh=dfdrug[(dfdrug['Compound']==compound) | (dfdrug['Compound'].str.contains('Vehicle'))]
                    
                    
                    
                    dfdrugLi=dfdrugVehi[(dfdrugVehi['Compound']==compound) | 
                                        (dfdrugVehi['Compound'].str.contains('Vehicle'))]
                    
                    

                    resdf=dfdrugVeh.groupby(['Experiment'])
                    resdfcomp=dfcompound.groupby(['Experiment'])
                    
                    
                    
                                          

                    if 'ISI' in df.columns:
                        tosub=folderISI

                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue
                                                                            (g, 'Compound', 'Group', intra=True, series=False), g, 
                                                        'Compound', 'ISI', 'Group', to, tosub, 'box', 'Label', hue_order_all, ageplot)))

                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue
                                                                                (g, 'Compound', 'Label', intra=True), g, 'Compound', 'ISI', 'Label', to, tosub, 'box', 'Group', hue_order_all, ageplot)))
                    if Parameter in df.columns:

                        tosub=folderMFR
                        #print('MFRplotting')
                        if (Parameter=='MFR(Hz)') or (Parameter=='BFR(Hz)'):
                            
                            if Li==True:
                                
                               
                            
                           
                                figcell, axcell=plt.subplots(figsize=(10, 8))
                                
                                sns.relplot(data=dfdrugLi,
                                        x="Age", y=Parameter,
                                        hue="Compound", hue_order=['Drugs', 'Vehicle'], style="Compound", style_order=['Drugs', 'Vehicle'], col="Group",
                                        height=4, aspect=.7, kind="line", ax=axcell)
                                plt.ylim([0, 10])
                               
                                    
                                plt.savefig(tosub+'/'+'Liplots'+'.png', dpi=300, bbox_inches='tight')
                                
                                

                                
                            #print(hue_order_all, 'hue_orders???')
                        
                        
                            fg=resdf.apply(lambda m : m.groupby(['LightWave']).apply(lambda g : 
                                                    relative(g, tosub,  hue_order_all, Parameter, 'Group', compound, to, 'True')))



                            gf=resdf.apply(lambda m : m.groupby(['LightWave']).apply(lambda g
                                            : relative(g, tosub,  hue_order_all, Parameter, 'Label', compound, to, 'True')))

                        
                  
                        
                        
                        
                       
                        res=resdf.apply(lambda m : m.groupby(['LightWave']).apply(lambda x:
                                                            x.groupby('Group').apply(lambda g: stat_plot_annot
                                                                (Xhue(g, 'Compound', 'Label', intra=True, 
                                                                                 series=True), g, 'Compound', 
                                                Parameter, 'Label', to, tosub, 'box', 'Group', hue_order_all, ageplot))))
                        
                        
                        
                        
                        res=resdf.apply(lambda m : m.groupby(['LightWave']).apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot
                                                                           (Xhue(g, 'Compound', 'Group',
                    intra=True, series=False), g, 'Compound', Parameter, 'Group', to, tosub, 
                                                                            'box', 'Label', hue_order_all, ageplot))))
                        
                        #print('printing?')
                        
                        res=resdfcomp.apply(lambda m : m.groupby(['LightWave']).apply(lambda x: x.groupby('Group').apply(lambda g: 
                                    stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True, series=True),
                                    g, 'Compound', Parameter, 'Label', to, tosub, 'box', 'Group', hue_order_all, True))))
                        
                        
                        
                        
                        
                    if 'awsr(min)' in df.columns:

                        tosub=folderawsr
                        #print(df[:3], 'awsr')
                        res=resdfcomp.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound', 'Label', intra=True), g, 'Compound', 'awsr(min)', 'Label', to, tosub, 'bar', 'Group', hue_order_all, ageplot)))
                        res=resdf.apply(lambda x: x.groupby('Group').apply(lambda g: stat_plot_annot
                                                                           (Xhue(g, 'Compound', 'Label', intra=True), g, 
                                                                            'Compound', 'awsr(min)', 'Label', to, tosub, 'bar', 'Group', hue_order_all, ageplot)))
                        res=resdf.apply(lambda x: x.groupby('Label').apply(lambda g: stat_plot_annot(Xhue(g, 'Compound',
                                                                                                          'Group', intra=True, series=False), g, 'Compound', 'awsr(min)', 'Group', to, tosub, 'bar', 'Label', hue_order_all, ageplot)))



    if "Cell Line" in effects:
        
        

        for df in dfs:
            
            #df=Filt_group(df, filt_dict)
            
            dfpharm=df

            df=df[df['Compound']=='No Group']
            
            if len(df)>0:

                resdf=df.groupby('Experiment')

                if 'ISI' in df.columns:
                    tosub=folderISI
                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True, series=False), g, 'Group', 'ISI', 'Well ID', to, tosub, 'box', 'Group', hue_order_all, ageplot))

                if Parameter in df.columns:
                    
                    #print(np.unique(df['Well ID'].values), 'wellfiltered?')
                    
                    

                    tosub=folderMFR
                    
                    
                    
                    pairs=Xhue(df, 'Age', 'Group', intra=True, series=False)
                    
                    #print(pairs, 'pairs, what')

                    pairs=stat_plot_annot(pairs,
                            df, 'Age', Parameter, 'Group', to, tosub, 'box', 'Group', hue_order_all, False)

                    res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True, series=False), 
                                                              
                                                              g, 'Group', Parameter, 'Well ID', to, tosub, 'box', 
                                                              'Group', hue_order_all, ageplot))


                if 'awsr(min)' in df.columns:
                    
                    #print(df[:3], 'awsr')

                    tosub=folderawsr
                    #res=resdf.apply(lambda g: stat_plot_annot(Xhue(g, 'Group', 'Well ID', intra=True),
                                                              #g, 'Group', 'awsr(min)', 'Well ID', to, tosub, 'box', 'Group', hue_order_all, ageplot))




        
    return 'Done'


# In[ ]:





# In[94]:


def PSTH_Plotting(df, targetlight, savep, to, hue, hueorder):
    
    exp=df['Experiment'].values[0]
    grp=df['Group'].values[0]
    dlabel=df['Label'].values[0]
    comp=df['Compound'].values[0]
    
    df['Label']=df['Label'].apply(lambda x : 'Baseline' if x=='Control' else x)
    
    
    
    ##df is dataset containing all parameter features for each channel for each well
    ##PSTH s are for 50 ms windows for 10ms windows, first window is 40 window
    
    ##Validation holds if stimlation with targetlight and controllight light occured. 
    
    

   ##'Stim' column provides info on stimulation , on or off psthcols=df.columns[df.columns.str.contains('PSTH')].values
    psthcols=df.columns[df.columns.str.contains('PSTH')].values
    #print(psthcols)
    #psthmax=max([int(txt.split('H')[1]) for txt in psthcols if txt.split('H')[1].isdigit()==True ])
    psthcolsf=[txt for txt in psthcols if txt.split('H')[1][::-1].split('-')[0].split('.')[1].isdigit()==True]
    print( psthcolsf)
    psthcols800=[txt for txt in psthcolsf if float(txt.split('H')[1])<500]
    ticklabs=[float(txt.split('H')[1]) for txt in psthcolsf if float(txt.split('H')[1])<500]
    psthmax=max([float(txt.split('H')[1]) for 
                 txt in psthcols800 if txt.split('H')[1][::-1].split('-')[0].split('.')[1].isdigit()==True])
    
    cols=[col for col in df.columns if col not in psthcols800]
    print(psthcols800, 'pdths')
   
    
    
    
    #for cl in psthcolsf:
        
        

        #df[cl]=(df[cl].values+0.001)/(df['PSTH0'].values+0.001)
        
        
    
    uvals=np.unique(df[hue].values)
    
    n_colors = len(uvals)
    colours = ['green', 'blue', 'grey', 'orange'][::-1]#cm.rainbow(np.linspace(0, 1, n_colors))
    
    handlelist=np.zeros(len(uvals), dtype='object')
    
    pvaldict={}
    pvaldrug={}

    
    
    fig, ax=plt.subplots(figsize=(8, 6))
    
    sns.set(font_scale=1.2)
    
    collect_val=pd.DataFrame()
    
    for index, uval in enumerate(uvals):
        
        
        pvaldict[uval]={}
        
        
            
        
        
        ##plot control light vs df light for each cell group/experiment condition
        t1=df[(df['LightWave']==targetlight) & (df['Stim']=='On') & (df[hue]==uval)]
        
        collect_val=pd.concat([collect_val, t1], axis=0)
        
        
        
        ind=np.where(hueorder==uval)[0]
        
        print(hueorder)
        
        print(uval, 'uval', ind)
        
        print(ind, hueorder[ind-1])
        
        if ind!=0:
            
            #print(ind, hueorder[ind-1], 'huet0')
            t0=df[(df['LightWave']==targetlight) & (df['Stim']=='On') & (df[hue]==hueorder[ind-1][0])]
            
            pvaldrug[uval]={}
        
        for pcol in psthcolsf:
            
            if ind!=0:
                
                try:
                
                    pvaldrug[uval][pcol]=wilcoxon(t1[pcol].values, t0[pcol].values)[1]
                except:
                    pass
                
                
            
            if pcol!='PSTH0':
                
                try:
                    
                
                    pvaldict[uval][pcol]=wilcoxon(t1['PSTH0'].values, t1[pcol].values)[1]
                    
                except:
                    pass#print('notsure')
                
                
        

        meltedt=pd.melt(t1, id_vars=cols, var_name='PSTH', value_name='Peak')
        
        




        sns.pointplot(data=meltedt, x='PSTH', y='Peak', join=True, color=colours[index], label=uval, ax=ax)  #colours[index]


        #maxpsthwindow
        
        print(index, ind, uval, 'ind, uval')



        targetlight_patch=mpatches.Patch(color=colours[index], label=uval)

        handlelist[ind[0]]=targetlight_patch



    ax.set_title(grp+dlabel)
    ax.tick_params(axis='x', labelrotation=90)

    ticklabels=ticklabs#np.arange(-25, psthmax-25, 5)
    plt.xticks(ticks=ax.get_xticks()[::5], labels=ticklabels[::5], fontsize=15)
    #plt.xlim([7, 90])
    #ax.set_xticks(ticklabels)
    
    #print(ticklabels, ax.get_xticks())
    #ax.set_yticks(fontsize=16)

    xspan=ax.get_xticks()

    ax.set_xlabel('PSTH')


    ax.set_ylabel('Relative Counts')
    
    plt.ylim([-1, 25])

    #ax.legend(handles=handlelist.tolist(), bbox_to_anchor=(1, 1))
    #plt.xlim([-100, 800])

    plt.savefig(to+'\\'+savep+exp+targetlight+grp+comp+dlabel+hue+'PSTHrbaseshortwowell5'+'.png')

    plt.show()
    
    pd.DataFrame(pvaldict).to_csv(to+'\\'+savep+exp+grp+comp+dlabel+hue+'StatPSTHtobase.csv')
    
    pd.DataFrame(pvaldrug).to_csv(to+'\\'+savep+exp+grp+comp+dlabel+hue+'StatPSTHtoPhase.csv')
    
    
    return collect_val


# add stim section and plot tiny box plots for psths

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[95]:


def age(ID, Unit):
    
    age='0'
    
   
    
    
   
    if Unit=='DIV':
    
        ageindex=ID.find(Unit)
        
        for d in np.arange(ageindex, len(ID), 1):
            if ID[d].isdigit():
                age=age+ID[d]
    else: 
        ageindex=ID.find(Unit)
        for d in np.arange(ageindex, len(ID), 1):
            if ID[d].isdigit():
                age=age+ID[d]
        
    
    return str(int(age[:]))
            

    
    


# In[96]:


def adjust_box_widths(g, fac):
    """
    Adjust the widths of a seaborn-generated boxplot.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])


# In[97]:


def segment_1(series, stamp, Label, Type, Duration=None): #??????
    
    if Type=='Spike':
        
        start=series.loc[series['Dose Label']==Label]['Timestamp [µs]'].min()
        end=series.loc[series['Dose Label']==Label]['Timestamp [µs]'].max()
        #print('end', end)
        #print('start', start )
        
    elif Type=='Burst':
        
        start=series.loc[series['Dose Label']==Label]['Timestamp [µs]'].min()
        end=series.loc[series['Dose Label']==Label]['Timestamp [µs]'].max()
    if Duration:
        
        end=start+(Duration*1000000)
    
    return (start, end)
        


# In[98]:


def segment(series, stamp, Label, Type, Duration=None):
    
    #print(Label, 'Label')#??????
    
    if Type=='Spike':
        
        start=int(stamp[str([Label])]['start'])
        end=int(stamp[str([Label])]['stop'])
        #print('end', end)
        #print('start', start )
        
    elif Type=='Burst':
        start=int(stamp[str([Label])]['start'])
        end=int(stamp[str([Label])]['stop'])
        #print('end', end)
        #print('start', start )
        
        
    if Duration:
        
        end=start+(Duration*1000000)
    
    return (start, end)
        


# In[99]:


#there can not be default params if list passed as an argument to a function, =X syntax is false
#neeed to really pass the variables


# In[100]:


def analyze(folder_path, identifier1, identifier2, to, winsize, Type, group_iN,
            effects, user_pharm, experimenttoassign, fontscale, ylabelsize,
            xlabelsize, figuresizex, figuresizey, overview, filt_dict, Unit='DIV', 
            equal=True, ISI_well=False, MC=True, Filter_Channel=True, 
            plot=True, assign=False, upper=False, acf=False, filt=True, scatter=False):
    
    
   
    
    #metadata files
    
    WellID=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellID.csv.npy")
    
    Welllabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellLabel.csv.npy")
    Channellabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\ChannelLabel.csv.npy")
    stamp_dict=np.load(folder_path+'/stamp.npy', allow_pickle=True)[()]
    led_dict=np.load(folder_path+'/LED_info.npy', allow_pickle=True)[()]
    
    #Data
    
    if Type=='Spike':
        f=load_exported(folder_path, identifier1, identifier2)[0]
    else:
        f=load_exported(folder_path, identifier1, identifier2)[1]
        
    #plotting configuration 
    plt.rcParams.update({'figure.max_open_warning': 0})
    
    sns.set(font_scale = fontscale)
    
    plt.rcParams['ytick.labelsize'] = ylabelsize
    plt.rcParams['xtick.labelsize'] = xlabelsize
    plt.rcParams["figure.figsize"] = (figuresizex, figuresizey)
    
    plt.rcParams['font.size']=16
    
    
    
    ##experiment name position 
    
    
    #Burst analysis
    
    
    
    #Variables, Metadata
    assignements=[] #placeholder for the group assignements 
    IDs=np.unique(f['Experiment'].values) #experiments 
    #f['Dose Label']=f['Dose Label'].apply(lambda x: '0' if x=='Control' else str(x)) #?? #dose label in str? already string?
    
    
    df_m=pd.DataFrame({'MFR(Hz)':[], 'Experiment':[], 'Compound':[], 'Stim':[], 'Label':[], 'Well ID':[], 'Group':[], 'Rchange':[],
                        'Channel ID':[]})
    
    df_ARSR=pd.DataFrame(columns=['awsr(min)', 'Well ID', 'Group', 'z score', 'Experiment', 'Compound', 'Label', 'Light', 'StimWells'])
    
    df_ISI=pd.DataFrame(columns=['ISI', 'Experiment', 'Compound', 'Label', 'Well ID', 'Group', 'Stim', 'StimWells'])
    
    hue_order_all={} 
    
    
         
    for ID in IDs:
        
        
        drugs_isi={}
        
        fileID=ID
        compounds=np.unique(f.loc[f['Experiment']==fileID]['Compound ID'].values)
        stamp_files=np.unique(f[f['Experiment']==fileID]['File'])
        stamp_filesorder=np.unique(f[f['Experiment']==fileID]['File'])
        stamp={}
       
     
        stamp[ID]={}
       
        for stamp_file in stamp_files:
            
            positionh5=stamp_file.index('h5')+2
           
            stamp.update(stamp_dict[stamp_file[:positionh5]])
      
        
        #print(stamp, 'new_stamp_file')
        #stamp=stamp_dict[stamp_file[:positionh5]]
        ##print(stamp, 'stamp')
        hue_order={}
        #as if all files for given experiment have same stamps
        
        #print('compound', compounds)
        
        
        ch_max=f['Channel ID'].value_counts().max()
        
        scale=ch_max*40000/120000000
        
    
        
       
        #fig1, axs1 = plt.subplots(24, len(compounds), figsize=(100, 15))
        
        for ind, compound in enumerate(compounds):
            
            fig=plt.figure(figsize=(16, 12))
            
            ax = fig.add_subplot(1,1,1)
            ax.set_title('Binned Spiking activity, drugs')
            ax.set_xlabel('Bins')
            ax.set_ylabel('N of spikes')
            
            
            
            figbin, axbins=plt.subplots(figsize=(10, 8))
            
            axbins.set_xlabel('Time (microsec)', fontsize=20)
            axbins.set_ylabel('Spike Count', fontsize=20)
            axbins.tick_params(axis='x', labelsize=20) 
            axbins.tick_params(axis='y', labelsize=20)
            
            
            
           
            
            
            #print(compound, 'next comp')
            
            
            
            
            
            hue=['Control']
            
            
            
            single_exp={'MFR(Hz)':[],  'Experiment':[], 'Compound':[], 'Stim':[], 'Label':[], 
                        'Well ID':[], 'Group':[], 'Rchange':[], 'Channel ID':[]}
            labels=np.unique(f[(f['Experiment']==fileID) & (f['Compound ID']==compound)]['Dose Label'].sort_values(ascending=True).values)
            series=f[(f['Experiment']==fileID) & (f['Compound ID']==compound)]
            led_condition_control=f[(f['Experiment']==fileID) & (f['Compound ID']==compound) & 
                                    (f['Dose Label']=='Control')]['File'].values[0]
            
            positionh5=led_condition_control.index('h5')+2
            
            
            wellsplot=np.unique(series['Well ID'].values)
            
            if (len([i for i in list(led_dict.keys()) if led_condition_control[:positionh5] in i]))>0:
                
                if len(list(led_dict[led_condition_control[:positionh5]].keys()))>0:
                    
                    ld=True
                else:
                    ld=False
            else:
                ld=False
            
            
            ###burst analysis
            if Type=='Burst':
        
                bz=Burst_analysis(to, f, group_iN, ld, l_amp)
            
            if overview=='True':
                if ld==True:
                    
                    plots=[plt.subplots(4, len(labels), figsize=(len(labels)*10, 50),
                                        sharey='row', gridspec_kw={'height_ratios': [2, 1, 1, 1]})
                           for  i in range(len(wellsplot))]
                    
                else:
                    plots=[plt.subplots(3, len(labels), figsize=(len(labels)*10, 30),
                                        sharey='row', gridspec_kw={'height_ratios': [2, 1, 1]})
                           for  i in range(len(wellsplot))]
                    
                    
                
                
            else:
                plots=[]
   
            labels=labels.astype('str')
            
            colors =['C{}'.format(i) for i in range(len(labels))]
            
            d_color=dict(zip(labels, colors))
            
            labord={}
            for label in labels:
                seg2=segment(series, stamp, label, Type)
                
                labord[label]=seg2[0]
            
            
            labord=sorted(labord.items(), key=lambda x:x[1])
            labord=np.array([i[0] for i in labord])
            #print(labord, 'labord')
            
                
            
            
            
            #lab_dict={'Control':0, 'AMPA1':2, 'CPP10':3, 'NBQX30':4, 'Vehicle':1}
            
            
          
           
          
            
            #print('label', labels)
            
            seg=segment(series, stamp, 'Control', Type)
            seg2=segment(series, stamp, 'Control', Type)
            #find start and end  of the label
            start, end=[seg[0], seg[1]]
            s, e=[seg2[0], seg2[1]]
            order=[s]
            dur=(end-start)/1000000
            
            
            
            
            
            if ld==True: #different function stimulated and baseline activity
                
                Led=led_dict[led_condition_control[:positionh5]]
                
                stat=mfr_bfr_Stim(to, fileID, compound, series, 50, winsize, start, end, s, e, plots, 0, 
                                   ch_max,  'Control', df_ISI, Type, acf, group_iN, Led, False)
                
                
                                  
                
                w=stat[3].astype(int)
                control=stat[1][w, :].ravel()
                control2=stat[9][w, :].ravel()
                single_exp=extension(single_exp, stat, control, control2, ID, compound,group_iN, 'Control')
                LED_list=stat[7]
                LED_IDs=stat[8]
                SelectedWells=[]
                ARSRresult=ARSR_apply(ID, to, group_iN, compound, 'Control', series, ld, LED_list, LED_IDs, df_ARSR, start, end, [], filt)
                SelectedWells=ARSRresult[1]
                df_ARSR=ARSRresult[0]
                #control_bins=stat[5]
                #axbins=plot_bins(control_bins, axbins, compound, 'Control', d_color['Control'])
                if ISI_well==True:
                    df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)
            else: 
                
                stat=mfr_bfr_label(to, fileID, compound, series, 50, winsize, start, end, s, e, plots, 0, 
                                   ch_max,  'Control', df_ISI, Type, acf, group_iN, scale, False)
                                  
               
                w=stat[3].astype(int) 
                control=stat[1][w, :].ravel()
                control2=stat[1][w, :].ravel()
                single_exp=extension(single_exp, stat, control, control2, ID, compound, group_iN, 'Control')
                LED_IDs=[]
                LED_list=[]
                ARSRresult=ARSR_apply(ID, to, group_iN, compound, 'Control', series, ld, LED_list, LED_IDs, df_ARSR, start, end, [], filt)
                SelectedWells=ARSRresult[1]
                df_ARSR=ARSRresult[0]
                #control_bins=stat[5]
                #axbins=plot_bins(control_bins, axbins, compound, 'Control', d_color['Control'])
                if ISI_well==True:
                    df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)
            
            plots=stat[6]
            #size=np.mean(stat[5], axis=0).size
            
            #ax.plot(np.arange(0, size, 1), np.mean(stat[5], axis=0), color=d_color['Control'], label='Control') 
            
            drugs_isi['Control']=stat[0]
            
            
            for label in labels:
                
             
                
                
                 if ((label !='0') and (label !='Control')):
                
                    seg=segment(series, stamp, label, Type)
                    seg2=segment(series, stamp, label, Type)
                    
                    start, end=[seg[0], seg[1]]
                    order.append(seg2[0])
                    hue.append(label)
                    dur=(end-start)/1000000
                    
                    led_condition=series[series['Dose Label']==label]['File'].values[0]
                    
                    positionh5=led_condition.index('h5')+2
                    if (len([i for i in list(led_dict.keys()) if led_condition[:positionh5] in i]))>0:
                        if len(list(led_dict[led_condition[:positionh5]].keys()))>0:
                            ld=True
                        else:
                            ld=False

                    else:
                        ld=False


                
               
                    
                    if ld==True:
                        Led=led_dict[led_condition[:positionh5]]

                        stat=mfr_bfr_Stim(to, fileID, compound, series, 50, winsize, start, end, s, e, 
                                           plots,  np.where(labord==label)[0][0], 
                                   ch_max,  label, df_ISI, Type, acf, group_iN, Led, False)
                        w=stat[3].astype(int)
                        single_exp=extension(single_exp, stat, control, control2, ID,compound, group_iN, label)
                        LED_list=stat[7]
                        LED_IDs=stat[8]
                        df_ARSR=ARSR_apply(ID, to, group_iN, compound,label, series, ld, LED_list, LED_IDs, df_ARSR, start, end, SelectedWells, filt)[0]
                        #label_bins=stat[5][:, :, :control_bins.shape[-1]]
                        #axbins=plot_bins(label_bins, axbins, compound, label, d_color[label])
                        if ISI_well==True:
                            df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)
                    
                        


                    else:

                        stat=mfr_bfr_label(to, fileID, compound, series, 50, winsize, start, end, s, e, 
                                           plots, np.where(labord==label)[0][0], ch_max, label, df_ISI, Type, acf, 
                                           group_iN, scale, False) 
                        
                
                        w=stat[3].astype(int)
                        single_exp=extension(single_exp, stat, control, control2, ID, compound,group_iN, label)
                        LED_list=[]
                        LED_IDs=[]
                        df_ARSR=ARSR_apply(ID, to, group_iN, compound, label, series, ld, LED_list, LED_IDs, df_ARSR, start, end, SelectedWells, filt)[0]
                        #label_bins=stat[5][:, :, :control_bins.shape[-1]]
                        #axbins=plot_bins(label_bins, axbins, compound, label, d_color[label])
                        if ISI_well==True:
                            df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)
                            
                    #si=np.mean(stat[5], axis=0).size
                    #if si>size:
                        
                        
                        #plt.plot(np.mean(stat[5], axis=0)[:size], color=d_color[label], label=label, figure=fig)
                    
                   
                    #else:
                        #pad_to=np.zeros(size-si).tolist()
                        
                        #ax.plot(np.arange(0, size, 1), (np.mean(stat[5], axis=0).tolist()+pad_to), color=d_color[label], label=label)
                        
                    plots=stat[6]  
                    
                    drugs_isi[label]=stat[0]
                    
                    #for row in range(stat[5].shape[0]):
                        
            axbins.legend()            #axs1[row, ind].plot(stat[5][row, :])
                    
            figbin.savefig(to+'/'+ID+compound+'BinsGroup'+'.png')                
            hue_order[compound]=sorted(list(zip(hue, order)), key = lambda x: x[1])
            
                
            df=pd.DataFrame(single_exp)
            
            if scatter==True:
            
                pcatter=scatter_doses(pd.DataFrame(single_exp), ID, to, np.array(labord))
            
            
           
            if ld==True:
                
                ledx1=int(round(LED_list[0]/winsize))-int(round(start/winsize))
                ledx2=int(round(LED_list[1]/winsize))-int(round(start/winsize))
                xbins=np.arange(ledx1, ledx2, 1)
                
                axbins.axhline(y=1.5*np.amax(stat[5]), xmin=ledx1, xmax=ledx2, color='blue')
                
                            
                            
            
              
                
            
            
            sj=[plots[i][0].savefig(to+'/'+str(wellsplot[i])+' '+compound+ID+'Well Overview'+'.tiff') for i in range(len(plots))]
           
       
            
           
            df=channel_select(df, 'MFR(Hz)', 'Channel ID', 'Group', MC, Filter_Channel, upper)
          
            df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)
            
       
            
        #fig1.savefig(to+'/'+ID+'BinsWells'+'png')
        
        hue_order_all[ID]=hue_order
            
    for df in [df_m, df_ARSR, df_ISI]:
        
        #print(df[:3], 'experiment')
        
        #adding age variable to datasets, constructing from ID
        
        df['Age']=df['Experiment'].apply(lambda x : age(x, Unit))
            
    tofile=to+'/'+'Files'  #Folder in to save data files

            
    if not os.path.exists(tofile):
        
        os.makedirs(tofile)   
    #saving data files to the folder    
    #df_ISI.to_excel(tofile+'/'+'ISI'+'.xlsx')
    df_m.to_excel(tofile+'/'+'MFR'+'.xlsx', )
    df_ARSR.to_excel(tofile+'/'+'AWSR'+'.xlsx')
    
           
  
           
    Res_string='Done unplotted' #can stop here
    
    #otherwise user defined part and effects to plot (but for r awsr, ISI, MFR), #f
   
    
    if plot==True:
        
        Res_string=Plot([df_m, df_ARSR, df_ISI], to, filt_dict, effects, hue_order_all)
        
    if assign==True:
        
        assignements=group_assign_II(df_ARSR[df_ARSR['Experiment']==experimenttoassign], user_pharm, equal)
    #savefilesintofilefolder
    #print(Res_string)
           
    return Res_string, assignements, stat[5]
                            
                            
        
                    
        


# In[101]:


def Meta(folder_path, identifier1, identifier2, to, winsize, Type, group_iN, Comps,
            effects, user_pharm, experimenttoassign, fontscale, ylabelsize,
            xlabelsize, figuresizex, figuresizey, overview, filt_dict, laborder, Unit='DIV', 
            equal=True, ISI_well=False, MC=True, Filter_Channel=True, 
            plot=True, assign=False, upper=False, acf=False, filt=True, scatter=False, aws=False, ageplot=True, match=False,
         newpharm=False, Li=True, osuc=True):
    
    
   
    
    #metadata files
    
    l_amp=[4, 6, 10, 12, 14, 15, 16, 18, 20]
    
    WellID=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellID.csv.npy")
    
    Welllabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellLabel.csv.npy")
    Channellabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\ChannelLabel.csv.npy")
    #if Type=='Burst':
        #inddex=folder_path.find('ursts')
        #mdata=folder_path[:inddex-1]+'csv'
        
    #else:
    mdata=folder_path
        
    stamp_dict=np.load(mdata+'/stamp.npy', allow_pickle=True)[()]
    led_dict=np.load(mdata+'/LED_info.npy', allow_pickle=True)[()]
    
    #Data
    
    if Type=='Spike':
        Params=['MFR(Hz)']
        
        #print(newpharm, 'newpharm')
        
        f=load_exported(folder_path, identifier1, identifier2, newpharm, Comps)[0]
    else:
        Params=['BFR(Hz)', 'Spike Count', 'Duration', 'Min ISI', 'Burstiness', 'Surprise']
        
        f=load_exported(folder_path, identifier1, identifier2, newpharm, Comps)[1]
        
        
    
        
    #plotting configuration 
    plt.rcParams.update({'figure.max_open_warning': 0})
    
    sns.set(font_scale = fontscale)
    
    plt.rcParams['ytick.labelsize'] = ylabelsize
    plt.rcParams['xtick.labelsize'] = xlabelsize
    plt.rcParams["figure.figsize"] = (figuresizex, figuresizey)
    
    plt.rcParams['font.size']=16
    
    
    
    ##experiment name position 
    
    
    #Burst analysis
    
    
    
    #Variables, Metadata
    assignements=[] #placeholder for the group assignements 
    IDs=np.unique(f['Experiment'].values) 
    Files=np.unique(f['File'].values)
    #experiments 
    #f['Dose Label']=f['Dose Label'].apply(lambda x: '0' if x=='Control' else str(x)) #?? #dose label in str? already string?
    
    #print(IDs, 'ids')
    
    
    
    iter_param=1
    
    for Param in Params:
        
        #print('Param', Param)
        
        
        if iter_param==1:
            
            if (Param=='MFR(Hz)') or (Param=='BFR(Hz)'):

                function='mfr-bfr'
                
                
                df_m=pd.DataFrame({Param:[], 'Experiment':[], 'Compound':[], 'Stim':[],
                                   'Label':[], 'Well ID':[], 'Group':[], 'Rchange':[],
                        'Channel ID':[], 
                                  'STTCpre':[],'LightWave':[]})
                
       
                df_ARSR=pd.DataFrame(columns=['awsr(min)', 'Well ID', 'Group', 'z score', 'Experiment', 
                                              'Compound', 'Label', 'Light', 'StimWells'])

                df_ISI=pd.DataFrame(columns=['ISI', 'Experiment', 'Compound', 'Label', 'Well ID', 'Group', 
                                             'Stim', 'StimWells'])

                hue_order_all={}
                
            else:
                
                
                function='dist-param'
                
                Param=[p for p in Params if p!='BFR(Hz)'] ###Param is here list of params
                
                dict_df={'Channel ID':[], 'Well ID':[], 'Well Label':[], 'Channel Label':[], 'Label':[],
       'Duration':[], 'Spike Count':[], 'Max ISI':[], 'Min ISI':[], 'Mean ISI':[], 'Variance':[],
       'Burstiness':[], 'Compound':[], 'Experiment':[], 'FT Spikes':[],
       'Start timestamp [µs]':[], 'Surprise':[], 'Stim':[], 'LightWave':[]}
                
                
                
                iter_param=0
                #print(Param, 'Parambursts')
                
                for p in Param:
                    
                    dict_df[p]=[]
                
                
                
                df_m=pd.DataFrame(dict_df)

               
            hue_order_all={}
                
                
            

            for ID in IDs:


                drugs_isi={}
                #ID=np.unique(f[f['File']==file]['Experiment']) #to be deleted

                fileID=ID
                compounds=np.unique(f.loc[f['Experiment']==fileID]['Compound ID'].values)
                doses=np.unique(f.loc[f['Experiment']==fileID]['Dose Label'].values)
                
               

                stamp_files=np.unique(f[f['Experiment']==fileID]['File'])
                stamp_filesorder=np.unique(f[f['Experiment']==fileID]['File'])
                stamp={}
                if newpharm==False:


                    stamp[ID]={}

                    for stamp_file in stamp_files:

                        positionh5=stamp_file.index('h5')+2

                        stamp.update(stamp_dict[stamp_file[:positionh5]])
                        #print(stamp, 'new_stamp_file')
             
                    


                
                #stamp=stamp_dict[stamp_file[:positionh5]]
                ##print(stamp, 'stamp')
                hue_order={}
                #as if all files for given experiment have same stamps

                #print('compound', compounds)


                ch_max=f['Channel ID'].value_counts().max()
                #print(ch_max, 'inmeta', f['Experiment'].value_counts())

                scale=ch_max*40000/120000000




                #fig1, axs1 = plt.subplots(24, len(compounds), figsize=(100, 15))
                
                if newpharm==True:
                    
                    
                    stamp={}

                    for dose in doses:

                        filename=f[(f['Experiment']==fileID) & (f['Dose Label']==dose)]['File'].values[0]

                       
                            
                        positionh5=filename.index('h5')+2
                        
                        

                        stamp[str([dose])]=stamp_dict[filename[:positionh5]][list(stamp_dict[filename[:positionh5]].keys())[0]]
                        
                        
                #print(stamp, 'bewstamp')
                        
                for ind, compound in enumerate(compounds):  
                    
                    
                    #print(compound, 'next comp')





                    hue=[]
                    
                    if function=='mfr-bfr':



                        single_exp={Param:[], 'Experiment':[], 'Compound':[], 'Stim':[],
                                   'Label':[], 'Well ID':[], 'Group':[], 'Rchange':[],
                        'Channel ID':[], 
                                  'STTCpre':[], 'LightWave':[]}
                
                        
                    else:
                        single_exp=dict_df
                        
                    labels=np.unique(f[(f['Experiment']==fileID) & (f['Compound ID']==compound)]
                                     ['Dose Label'].sort_values(ascending=True).values)
                    series=f[(f['Experiment']==fileID) & (f['Compound ID']==compound)]
                    
                    if 'Control' in labels:
                        
                        led_condition_controls=np.unique(f[(f['Experiment']==fileID) & (f['Compound ID']==compound) & 
                                                (f['Dose Label']=='Control')]['File'].values)
                    else:
                        led_condition_controls=np.unique(f[(f['Experiment']==fileID) & (f['Compound ID']==compound)]
                                                         ['File'].values)
                    
                    wellsplot=np.unique(series['Well ID'].values)
                    
                    
                    if ((overview=='True') and (function=='mfr-bfr')):
                        if len(led_condition_controls)>0:

                            plots=[plt.subplots(5, len(labels)+len(led_condition_controls)-1, figsize=
                                                ((len(labels)+len(led_condition_controls)-1)*25, 50),
                                                sharey='row', gridspec_kw={'height_ratios': [2, 1, 1, 1, 1]})
                                   for  i in range(len(wellsplot))]

                        else:
                            plots=[plt.subplots(3, len(labels), figsize=(len(labels)*25, 30),
                                                sharey='row', gridspec_kw={'height_ratios': [2, 1, 1]})
                                   for  i in range(len(wellsplot))]




                    else:
                        plots=[]
                        
                    lediter=0
                        
                    #print(led_condition_controls, 'looper')
                    for  led_condition_control in  led_condition_controls:

                        positionh5=led_condition_control.index('h5')+2




                        #print(len(wellsplot), 'wellsplot')

                        #print(led_condition_control[:positionh5], 'led_conditions')

                        if (len([i for i in list(led_dict.keys()) if led_condition_control[:positionh5] in i]))>0:
                            
                            if len(list(led_dict[led_condition_control[:positionh5]].keys()))>0:
                                
                                key0=list(led_dict[led_condition_control[:positionh5]].keys())[0]
                                
                                if len(led_dict[led_condition_control[:positionh5]][key0])>0:
                                    

                                    ld=True
                                else:
                                    ld=False
                            else:
                                ld=False
                        else:
                            ld=False

                        ###burst analysis




                        labels=labels.astype('str')

                        colors =['C{}'.format(i) for i in range(len(labels))]

                        d_color=dict(zip(labels, colors))

                        labord={}
                        for label in labels:
                            seg2=segment(series, stamp, label, Type)

                            labord[label]=seg2[0]


                        labord=sorted(labord.items(), key=lambda x:x[1])
                        labord=np.array([i[0] for i in labord])
                        #print(labord, 'labordfirst')





                        #lab_dict={'Control':0, 'AMPA1':2, 'CPP10':3, 'NBQX30':4, 'Vehicle':1}






                        #print('label', labels)

                       

                        if newpharm==True:

                            if len(laborder)>0:

                                labord=np.array(laborder)

                                #print(labord, 'labordii')
                                
                        
                        order=[]

                        if 'Control' in labels:
                            hue.append('Control')
                            seg=segment(series, stamp, 'Control', Type)
                            seg2=segment(series, stamp, 'Control', Type)
                            #find start and end  of the label
                            start, end=[seg[0], seg[1]]
                            s, e=[seg2[0], seg2[1]]
                            order=[s]
                            dur=(end-start)/1000000
                            
                            print(dur, 'duration', start, end, s, e)



                            if ld==True: #different function stimulated and baseline activity

                                Led=led_dict[led_condition_control[:positionh5]]

                                #print(led_condition_control[:positionh5], 'led_info')

                                if function=='mfr-bfr':

                                    stat=mfr_bfr_Stim(to, fileID, compound, series, 50, winsize, start, end, s, e, plots, 0+lediter, 
                                                       ch_max,  'Control', df_ISI, Type, Param, acf, group_iN, Led, 
                                                      led_condition_control[:positionh5], False, osuc,  l_amp)


                                    w=stat[3].astype(int)
                                    control=stat[1][w, :].ravel()
                                    control2=stat[9][w, :].ravel()
                                    single_exp=extension(single_exp, stat, control, control2, ID, compound,group_iN, 
                                                         'Control', Param, function, True)
                                    LED_list=stat[7]
                                    LED_IDs=stat[8]
                                    SelectedWells=[]
                                    if aws==True:
                                        ARSRresult=ARSR_apply(ID, to, group_iN, compound, 'Control', series, ld, 
                                                              LED_list, LED_IDs, df_ARSR, start, end, [], filt)
                                        SelectedWells=ARSRresult[1]
                                        df_ARSR=ARSRresult[0]
                                    #control_bins=stat[5]
                                    #axbins=plot_bins(control_bins, axbins, compound, 'Control', d_color['Control'])
                                    if ISI_well==True:
                                        df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)
                                    plots=stat[6]

                                    drugs_isi['Control']=stat[0]
                                else:
                                    stat=Burst_analysis(to, ID, compound, series,  start, end, s, e, 
                                                        0, ch_max, 'Control', Type,  group_iN, Param, Led, ld, l_amp)



                                    single_exp=stat[0]
                                    df=pd.DataFrame(single_exp)

                                    if len(df_m['Channel ID'])<1:
                                        df_m=df

                                    else:
                                        df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)
                                    
                                    #print(np.unique(single_exp['Label'].values), '???Conntrol')


                                   

                                    ###rchange

                                    LED_list=stat[1]
                                    LED_IDs=stat[2]


                            else: 

                                if function=='mfr-bfr':

                                    #print('came here')

                                    stat=mfr_bfr_label(to, fileID, compound, series, 50, winsize, start, end, s, e, plots, 0, 
                                                       ch_max,  'Control', df_ISI, Type, acf, group_iN, scale, False)





                                    w=stat[3].astype(int) 
                                    control=stat[1][w, :].ravel()
                                    control2=stat[1][w, :].ravel()
                                    single_exp=extension(single_exp, stat, control, control2, ID, compound,group_iN, 'Control', 
                                                         Param, function, False)
                                    LED_IDs=[]
                                    LED_list=[]
                                    if aws==True:
                                        ARSRresult=ARSR_apply(ID, to, group_iN, compound, 'Control', series, ld, LED_list, LED_IDs, df_ARSR, start, end, [], filt)
                                        SelectedWells=ARSRresult[1]
                                        df_ARSR=ARSRresult[0]
                                    #control_bins=stat[5]
                                    #axbins=plot_bins(control_bins, axbins, compound, 'Control', d_color['Control'])
                                    if ISI_well==True:
                                        df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)

                                    plots=stat[6]

                                    drugs_isi['Control']=stat[0]


                                else:
                                    #print('burst here')
                                    stat=Burst_analysis(to, ID, compound, series,  start, end, s, e, 
                                                       0, ch_max, 'Control', Type,  group_iN, Param, [], ld, l_amp)



                                    single_exp=stat[0]
                                    
                                    df=pd.DataFrame(single_exp)

                                    if len(df_m['Channel ID'])<1:
                                        df_m=df

                                    else:
                                        df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)

                                   

                                    ###rchange

                                    LED_list=stat[1]
                                    LED_IDs=stat[2]


                        if len(laborder)>0:

                            labord=np.array(laborder)

                            labels=laborder



                        for label in labels:




                             if ((label !='0') and (label !='Control')):

                                seg=segment(series, stamp, label, Type)
                                seg2=segment(series, stamp, label, Type)

                                start, end=[seg[0], seg[1]]
                                order.append(seg2[0])
                                s, e=[seg2[0], seg2[1]]
                                hue.append(label)
                                dur=(end-start)/1000000

                                led_condition=series[series['Dose Label']==label]['File'].values[0]

                                positionh5=led_condition.index('h5')+2
                                if (len([i for i in list(led_dict.keys()) if led_condition[:positionh5] in i]))>0:
                                    if len(list(led_dict[led_condition[:positionh5]].keys()))>0:
                                        ld=True
                                    else:
                                        ld=False

                                else:
                                    ld=False





                                if ld==True:
                                    Led=led_dict[led_condition[:positionh5]]

                                    #print(led_dict[led_condition[:positionh5]], 'led_dict_label')

                                    if function=='mfr-bfr':

                                        #print(labord, label, start, end, 'labordlabeltestingcpp')

                                        stat=mfr_bfr_Stim(to, fileID, compound, series, 50, winsize, start, end, s, e, 
                                                           plots,  np.where(labord==label)[0][0]+lediter, 
                                                   ch_max,  label, df_ISI, Type, Param, acf, group_iN, Led, 
                                                          led_condition[:positionh5], False, osuc,  l_amp)
                                        w=stat[3].astype('int64')
                                        single_exp=extension(single_exp, stat, control, control2, 
                                                             ID,compound, group_iN, label, Param, function, True)
                                        LED_list=stat[7]
                                        LED_IDs=stat[8]
                                        if aws==True:
                                            df_ARSR=ARSR_apply(ID, to, group_iN, compound,label, series, ld, LED_list, LED_IDs, df_ARSR, start, end, SelectedWells, filt)[0]
                                        #label_bins=stat[5][:, :, :control_bins.shape[-1]]
                                        #axbins=plot_bins(label_bins, axbins, compound, label, d_color[label])
                                        if ISI_well==True:
                                            df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)
                                        plots=stat[6]  

                                        drugs_isi[label]=stat[0]
                                    else:
                                        stat=Burst_analysis(to, ID, compound, series,  start, end, s, e, 
                                                    np.where(labord==label)[0][0], ch_max, label, Type,  group_iN, 
                                                            Param, Led, ld, l_amp)




                                        single_exp=stat[0]
                                        
                                        df=pd.DataFrame(single_exp)

                                        if len(df_m['Channel ID'])<1:
                                            df_m=df

                                        else:

                                            df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)


                                       

                                        ###rchange

                                        LED_list=stat[1]
                                        LED_IDs=stat[2]




                                else:

                                    if function=='mfr-bfr':

                                        stat=mfr_bfr_label(to, fileID, compound, series, 50, winsize, start, end, s, e, 
                                                           plots, np.where(labord==label)[0][0], ch_max, label, 
                                                           df_ISI, Type,
                                                           acf, 
                                                           group_iN, scale, False) 


                                        w=stat[3].astype('int64')
                                        
                                        if 'Control' not in labels:
                                            control=stat[1][w, :].ravel()
                                            control2=stat[1][w, :].ravel()
                                            
                                        single_exp=extension(single_exp, stat, control, control2, ID,
                                                             compound,group_iN, label, Param, function, False)
                                        LED_list=[]
                                        LED_IDs=[]
                                        if aws==True:

                                            df_ARSR=ARSR_apply(ID, to, group_iN, compound, label, series, ld, LED_list, LED_IDs, df_ARSR, start, end, SelectedWells, filt)[0]
                                        #label_bins=stat[5][:, :, :control_bins.shape[-1]]
                                        #axbins=plot_bins(label_bins, axbins, compound, label, d_color[label])
                                        if ISI_well==True:
                                            df_ISI=pd.concat([df_ISI, stat[2]], axis=0).reset_index(drop=True)


                                        plots=stat[6]  

                                        drugs_isi[label]=stat[0]


                                    else:
                                        stat=Burst_analysis(to, ID, compound, series,  start, end, s, e, 
                                                    np.where(labord==label)[0][0], ch_max, label, Type,  group_iN,
                                                            Param, [], ld, l_amp)



                                        single_exp=stat[0]
                                        
                                        df=pd.DataFrame(single_exp)

                                        if len(df_m['Channel ID'])<1:
                                            
                                            df_m=df

                                        else:

                                            df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)


                                        

                                        ###rchange

                                        LED_list=stat[1]
                                        LED_IDs=stat[2]


                        lediter=lediter+1
                       


                        hue_order[compound]=sorted(list(zip(hue, order)), key = lambda x: x[1])
                        
                        if scatter==True:
                            
                            #improve for bursts as well
                            
                            
                            pcatter=scatter_doses(pd.DataFrame(single_exp), ID, to, np.array(labord))


                    if function=='mfr-bfr':
                        for key in list(single_exp.keys()):
                            
                            print(key, len(single_exp[key]))
                        
                        
                        
                        df=pd.DataFrame(single_exp)

                        if len(df_m['Channel ID'])<1:
                            
                            df_m=df

                        else:

                            df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)

                    hue_order_all[ID]=hue_order

                       
                    if function=='mfr-bfr':
                        
                        folderwells=to+'/'+'WellOvervew'
                        
                        if not os.path.exists(folderwells):
                            os.makedirs(folderwells)






                        sj=[plots[i][0].savefig(folderwells+'/'+str(wellsplot[i])+' '+compound+ID+'Well Overview'+'.tiff') for i in range(len(plots))]




                    #df=channel_select(df, Param, 'Channel ID', 'Group', MC, Filter_Channel, upper)
                    
                    

                        #df_m=pd.concat([df_m, df], axis=0).reset_index(drop=True)
                    
                    ##print(df_m.describe(), 'concate?')



                    #fig1.savefig(to+'/'+ID+'BinsWells'+'png')

                    
                    
            #print(hue_order_all, 'hueorder')
                
            if function=='mfr-bfr':

                for df in [df_m, df_ARSR, df_ISI]:

                    #print(np.unique(df['Experiment'].values), 'experiment')

                    #adding age variable to datasets, constructing from ID

                    df['Age']=df['Experiment'].apply(lambda x : age(x, Unit)  if x!='iN7_0_BP01.3BP29.4DIV6' else str(64))

                tofile=to+'/'+'Files'  #Folder in to save data files
                
            else:
                for df in [df_m]:
                    #print(df[:3], 'experiment')

                    #adding age variable to datasets, constructing from ID

                    df['Age']=df['Experiment'].apply(lambda x : age(x, Unit) if x!='iN7_0_BP01.3BP29.4DIV6' else str(64))
                tofile=to+'/'+'Files'  #Folder in to save data files
                


            if not os.path.exists(tofile):

                os.makedirs(tofile)   
            #saving data files to the folder    
            #df_ISI.to_excel(tofile+'/'+'ISI'+'.xlsx')
            
            if function=='mfr-bfr':
                df_m.to_excel(tofile+'/'+Param+'.xlsx')

                if aws==True:
                    df_ARSR.to_excel(tofile+'/'+'AWSR'+'.xlsx')
                    
                    
            else:
                #df_m=df_m.dropna() ###here produces uniqeal, 
                df_m.to_excel(tofile+'/'+'Burst_analysis'+'.xlsx')
                




            Res_string='Done unplotted' #can stop here

            #otherwise user defined part and effects to plot (but for r awsr, ISI, MFR), #f
            
            if match==True:

                #df_m=df_m.groupby('Group').apply(lambda x: match_analysis(x, Param)).reset_index(drop=True)
                df_m=df_m.groupby('Group').apply(lambda x: match_f(x, Param)).reset_index(drop=True)
            
            if plot==True:
                
                if function=='mfr-bfr':
                    
                    #print('here')
                    
                    #print(hue_order_all, 'hue_order_all')
                    np.save(to+'/'+'hue_order', hue_order_all)
                    
                    Paramses=[Param]
                    
                    for par in Paramses:
                        
                        #if par=='MFR(Hz)':
                            
                            #df_m=df_m[df_m['MFR(Hz)']>=0.1]
                        
                       

                        Res_string=Plot([df_m], to, filt_dict, effects, hue_order_all, ageplot, par, Li)
                        
                        
                else:
                    for p in Param:
                        
                      
                            
                            
                        
                        #print('BFR param plotting')
                        Res_string=Plot([df_m], to, filt_dict, effects, hue_order_all, ageplot, p, Li)
                        

            if assign==True:

                assignements=group_assign_II(df_ARSR[df_ARSR['Experiment']==experimenttoassign], user_pharm, equal)
            #savefilesintofilefolder
            #print(Res_string)

    return (Res_string, assignements, stat, plots) #stat[5]
                            
                            


# In[102]:


def hue_change(df, hue_order_all, compound):
    
    ide=df['Experiment'].values[0]
    
    hue_order=hue_order_all[ide][compound]
    hue_orderrel=[i[0] for i in hue_order]
    
    controlindex=np.where(np.array(hue_orderrel)=='Control')[0]
    
    hueplot=np.delete(np.array(hue_orderrel), int(controlindex))
    
    return (hue_orderrel, hueplot)
    
    


# In[103]:


def relative1(df, hue_order_all, Param,  groupby, compound, to):
    
    
    
  
    
    
    group=df['Group'].values[0]
    ide=df['Experiment'].values[0]
    
    #print(np.unique(df['Compound'].values), df.columns, 'fjksjfljshf')
    
    hue_order=hue_order_all[ide][compound]
    
    hue_order=[i[0] for i in hue_order]
    
    #print(hue_order, 'hueordrelative')
    
    
    
    
    
    labels=np.unique(df['Label'].values)
    
    dfnew=pd.DataFrame(columns=df.columns.tolist()+['Change'])
    
    for label in labels:
        
        if label!='Control':
            
            ind=np.where(np.array(hue_order)==label)[0]
            
            dfpre=df[df['Label']==hue_order[int(ind-1)]].sort_values(by='Channel ID')
            
            dfpost=df[df['Label']==hue_order[int(ind)]].sort_values(by='Channel ID')
            
            x=dfpre[Param].values
            y=dfpost[Param].values
            
            r=rchange(x, y)
            
            
            
            dfpost['Change']=r
            
            #print(dfpost[:3])
            
            dfnew=pd.concat([dfnew, dfpost], axis=0)
            
            
    
    pair=Xhue(dfnew, 'Compound', 'Label', intra=True, series=True)[1]        
    plt.figure(figsize=(10, 8))
    
    hue_orderplot=np.delete(np.array(hue_order), int(np.where(np.array(hue_order)=='Control')[0]))
    
    ##print(np.unique(dfnew['Compound'].values), 'newvjfn')
    
    ax = sns.boxplot(data=dfnew, x='Compound', y='Change', hue='Label', hue_order=hue_orderplot)#, order=order, hue_order=hue_order)

    sns.stripplot(data=dfnew, x='Compound', y='Change', hue='Label', hue_order=hue_orderplot, ax=ax, dodge=True)                 

    annot = Annotator(ax, pair, data=dfnew, x='Compound', y='Change',
                      hue='Label', hue_order=hue_orderplot)#, order=order, hue_order=hue_order)


    annot.configure(test='Kruskal', text_format='simple')
    annot.apply_test()
    annot.annotate()
    
    plt.ylabel('Change of'+' '+Param, fontsize=20)
    #plt.xlabel(fontsize=20)
    plt.yticks(ymin=0, fontsize=20)
    
    plt.savefig(to+'/'+ide+compound+group+'Change'+Param+'.tiff', bbox_inches='tight')

    #plt.legend('',frameon=False)

 




    
   


# In[104]:


def plot_radar(dt, compare, groupbie, Norm, rangesext, tosub):
    
    gp=np.unique(dt[compare].values)
    
    
    d=''
    for l in groupbie:
        d=d+np.unique(dt[l].values)[0]
        dt=dt.drop(l, axis=1)
   
    
    
    columns=dt.columns
    
    #print(columns)
    
    variables =dt.drop(compare,  axis=1).columns
    
    if Norm=='False':
        
        suf='Norm'
        
        ranges=rangesext
    else:
        
        suf='Sym'
        rmin=np.ones([len(variables)])*-1
        rmax=np.ones([len(variables)])
        ranges=list(zip(rmin, rmax))
        
        
    #print(len(variables), len(ranges))
    
    #print(dlabel)
    
    fig1 = plt.figure(figsize=(6, 6))
    
    #plt.title(dlabel)
                      
                      
    
    my_palette = plt.cm.get_cmap("PiYG", len(dt.index))
    
    handles=[]
    
    for ind, g in enumerate(gp): 
        data = dt[dt[compare]==g].drop(compare, axis=1).values.ravel()
        
        
        
        ##print(data, variables.values, ranges, my_palette(ind))
       
        radar = ComplexRadar(fig1, variables.values, ranges)
        radar.plot(data, g,  color=my_palette(ind), alpha=1)
        #radar.fill(data, color=my_palette(ind), alpha=0.4)
        r_patch = mpatches.Patch(color=my_palette(ind), label=g)
        handles.append(r_patch)
        
    plt.legend(handles=handles, loc='best', bbox_to_anchor=(0.1, 0))
    
   

    fig1.savefig(tosub+'/'+d+'Radarby'+compare+suf+'.svg', format='svg', bbox_inches='tight')
    plt.show()    
        


# In[105]:


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
        


# In[106]:


def MinMax(data, lim):
    """Lim is a featue range given by list [min, max]"""
    
    
    
    if len(np.unique(data))==1:
        
        
        data_scaled=data
    
    else:
        
        
        data_scaled=[(((lim[-1]-lim[0])*(x-np.min(data))/(np.max(data)-np.min(data)))+lim[0]) if np.unique(data).shape[0]!=1
                     else x
                     for x 
                     in data]
        
    
    
    return data_scaled


# In[107]:


class ComplexRadar():
    def __init__(self, fig, variables, ranges,
                 n_ordinate_levels=4):
        angles = np.arange(0, 360, 360./len(variables))

        axes = [fig.add_axes([0.1,0.1,0.9,0.9], polar=True,
                label = "axes{}".format(i))
    
                for i in range(len(variables))]
        
        
        l, text = axes[0].set_thetagrids(angles, 
                                         labels=[v[6:] if 'Change' in v else v for v in variables], fontsize=12)
        [txt.set_rotation(angle-90) for txt, angle 
             in zip(text, angles)]
        axes[0].tick_params(axis='both', which='major', pad=40)
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
        
        #print(data)
        sdata = _scale_data(data, self.ranges)
        
        ##print(self.angle, data, sdata[0])
        self.ax.plot(self.angle, np.r_[sdata, sdata[0]], label=g,  *args, **kw)
    def fill(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.fill(self.angle, np.r_[sdata, sdata[0]], *args, **kw)


# In[108]:


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


# In[ ]:





# In[109]:


def relative(df, tosub, hue_order_all, Param, groubby, compound, to, save):
    
    #print(df[:3], 'dgdg')
    
    if compound!='Control' and compound!='No Group':
        
    

        ide=df['Experiment'].values[0]
        #print(ide, 'ideee???')
        lightw=df['LightWave'].values[0]
        
        #print(hue_order_all[ide], 'hueorderide???')

       

        hue_order=hue_order_all[ide][compound]

        hue_order=[i[0] for i in hue_order]
        
        #print(hue_order, 'hue_order')

       

        labels=np.unique(df['Label'].values)

        
        dfnew=pd.DataFrame(columns=df.columns.tolist()+['Change'])

        for label in labels:

            if label!='Control':

                ind=np.where(np.array(hue_order)==label)[0]

                #print(ind, label, 'inlopdjgh', Param, 'Param')

                dfpre=df[df['Label']==hue_order[int(ind-1)]].sort_values(by='Channel ID')

                dfpost=df[df['Label']==hue_order[int(ind)]].sort_values(by='Channel ID')

                x=dfpre[Param].values
                y=dfpost[Param].values

                r=rchange(x, y)
                
                #print(len(r), 'length', len(dfpre), len(dfpost), label)



                dfpost['Change'+Param]=r

               

                dfnew=pd.concat([dfnew, dfpost], axis=0)
                
        if save=='True':
                
        

            dfnew.to_csv(to+'/Files'+'/Change of'+Param+'.csv') ## changed this ####
        if groubby=='Group':
            
            


            groups=np.unique(dfnew['Group'].values)

            #print('Groups', groups)

            for group in groups:

                #print('group', group)

                dfgp=dfnew[dfnew['Group']==group]


                if len(np.unique(dfgp['Label'].values))>1:
                    
                    ordlist=None

                    

                    #print(np.unique(dfgp['Compound'].values).tolist(), 'compound in relative')

                    if 'Vehicle' in np.unique(dfgp['Compound'].values).tolist():
                        order=['Vehicle']
                        ordlist=np.unique(dfgp['Compound'].values).tolist()
                        ordlist.remove('Vehicle')
                        #print(ordlist, 'ordlist')
                        ordlist.extend(order)

                        #print('ordercompound', order)
                    else:
                        order=None

                    #print(np.unique(dfgp['Compound'].values), 'uniquecompound')
                    
                    if len(np.unique(dfgp['Compound'].values).tolist())>1:
                        
                        pair=Xhue(dfgp, 'Compound', 'Label', intra=True, series=True)[1]
                        #print(pair, "????")
                        test1='Mann-Whitney'
                    else:
                        pair=Xhue(dfgp, 'Compound', 'Label', intra=True, series=True)[0]
                        
                        #print(pair, "????new")
                        test1='Wilcoxon'
                        
                        


                    plt.figure(figsize=(10, 8))





                    hue_orderplot=np.delete(np.array(hue_order), int(np.where(np.array(hue_order)=='Control')[0]))

                    #print(hue_order, 'hue_order', hue_orderplot, 'hue_orderplot')

                    ##print(np.unique(dfnew['Compound'].values), 'newvjfn')

                    #print(np.unique(dfgp[dfgp['Compound']=='Drugs']['Label'].values), 'lab')

                    ax = sns.boxplot(data=dfgp, x='Compound', y='Change'+Param, hue='Label', 
                                     hue_order=hue_orderplot, order=ordlist, showfliers=False)#, order=order, hue_order=hue_order)

                    sns.stripplot(data=dfgp, x='Compound', y='Change'+Param, hue='Label',
                                  hue_order=hue_orderplot, order=ordlist, ax=ax, dodge=True)                 

                    annot = Annotator(ax, pair, data=dfgp, x='Compound', y='Change'+Param,
                                      hue='Label', hue_order=hue_orderplot, order=ordlist)#, order=order, hue_order=hue_order)


                    annot.configure(test=test1, text_format='simple')
                    annot.apply_test()
                    annot.annotate()

                    plt.ylabel('Change of'+' '+Param, fontsize=20)
                    #plt.xlabel(fontsize=20)
                    #plt.yticks(ymin=0, fontsize=20)

                    plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))

                    plt.savefig(to+'/'+ide+str(lightw)+compound+group+'Change'+Param+'.tiff', bbox_inches='tight')
                    plt.show()

                #plt.legend('',frameon=False)



        else:

            Title='byLabel'
            labs=np.unique(dfnew['Label'].values)

            for lab in labs: 

                dfl=dfnew[dfnew['Label']==lab]
                
                if len(np.unique(dfl['Group'].values))>1:
                    pairs=Xhue(dfl, 'Compound', 'Group', intra=True, series=False)
                    #print(pairs, 'pairsrela')
                    pair=pairs[0]
                    ordlist=None
                    
                    if 'Vehicle' in np.unique(dfl['Compound'].values).tolist():
                        order=['Vehicle']
                        ordlist=np.unique(dfl['Compound'].values).tolist()
                        ordlist.remove('Vehicle')
                        #print(ordlist, 'ordlist')
                        ordlist.extend(order)
                        
                        #print('ordercompound', order)
                    else:
                        order=None
                        
                    #print(np.unique(dfl['Compound'].values), 'uniquecompound')

                        
                    #print('ordercom', order)

                    plt.figure(figsize=(10, 8))

                    hue_orderplot=np.delete(np.array(hue_order), int(np.where(np.array(hue_order)=='Control')[0]))

                    ##print(np.unique(dfnew['Compound'].values), 'newvjfn')

                    ax = sns.boxplot(data=dfl, x='Compound', y='Change'+Param, hue='Group', order=ordlist, showfliers=False)#, order=order, hue_order=hue_order)

                    sns.stripplot(data=dfl, x='Compound', y='Change'+Param, hue='Group', order=ordlist, ax=ax, dodge=True)                 

                    annot = Annotator(ax, pair, data=dfl, x='Compound', y='Change'+Param,
                                      hue='Group', order=ordlist)#, order=order, hue_order=hue_order)


                    annot.configure(test='Mann-Whitney', text_format='simple')
                    annot.apply_test()
                    annot.annotate()

                    plt.ylabel('Change of'+' '+Param, fontsize=20)
                    #plt.xlabel(fontsize=20)
                    #plt.yticks(ymin=0, fontsize=20)
                    plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1))

                    plt.savefig(to+'/'+ide+str(lightw)+compound+Title+lab+'Change'+Param+'.tiff', bbox_inches='tight')

                    plt.show()

                #plt.legend('',frameon=False)

        return dfnew

           

        





   


# In[110]:


def plot_bins(array, axbin, compound, label, color):
    
    array=array.reshape(array.shape[0]*array.shape[1], array.shape[2])
    
    
    
    array=np.array([row for row in array if len(row[row>0])>10])
    
    xscatter=np.arange(1, array.shape[-1]+1, 1)
    xscatter2d=np.repeat(xscatter, array.shape[0], axis=0)
    #print(xscatter.shape, 'shapex')
    #print(array.shape, 'arrayshape')
    
    
    #bp=axbin.boxplot(array.T.tolist(), showfliers=False, patch_artist=True, widths=0.25, showbox=False)
    
    
    #for patch in bp['boxes']:
        
        
        #patch.set(facecolor=color)
        #patch.set_alpha(0.5)
    
   
    axbin.plot(np.mean(array, axis=0), color=color, label=label)
    #axbin.scatter(xscatter2d, array, color=color)
    
        
        
    return axbin
    
   
    
    
    
    


# In[ ]:





# In[111]:


def ARSR_apply(ID, to, group_iN, compound, label, series, ld, LED_list, LED_IDs, df_ARSR, start, end, SelectedWells, filt):
    
    
    if ld==True:
        
        LED_s=LED_list[0]
        LED_e=LED_list[1]
        
    
        seriesawsr=series[series['Dose Label']==label]
        
        spikes_led=seriesawsr[(seriesawsr['Timestamp [µs]']>LED_s) & (seriesawsr['Timestamp [µs]']<LED_e)]
        
        spike_base=seriesawsr[seriesawsr['Timestamp [µs]']<LED_s]#baseline pre activity, in the future double condition not stimulated el and light

        spike_post=seriesawsr[seriesawsr['Timestamp [µs]']>LED_e]
            
        durled=(LED_e-LED_s)/1000000
        durbase=(LED_s-start)/1000000
        durpost=(end-LED_e)/1000000
        
        arsrbase=ARSR(spike_base, group_iN, durbase, filt).values[0]
        dfarsr=pd.DataFrame(arsrbase)
        dfarsr['awsr(min)']=dfarsr['awsr(min)'].apply(lambda x : np.log(x))

        #zscore well selection 

        dfarsr['z score']=np.abs(stats.zscore(dfarsr['awsr(min)']))

        #outlier well ids are defined in controls
        
        if len(SelectedWells)==0:

            SelectedWells=dfarsr[dfarsr['z score']<3]['Well ID'].values

        dfarsr=dfarsr[dfarsr['Well ID'].isin(SelectedWells)]

        dfarsr['Compound']=np.repeat(compound, len(dfarsr))
        dfarsr['Label']=np.repeat(label, len(dfarsr))
        dfarsr['Experiment']=np.repeat(ID, len(dfarsr))
        dfarsr['Light']=np.repeat('off', len(dfarsr))
        dfarsr['StimWells']=dfarsr['Well ID'].apply(lambda x: "on" if x in LED_IDs else 'off')
        
        
        df_ARSR=pd.concat([df_ARSR, dfarsr], axis=0).reset_index(drop=True)
        
        dfarsr=dfarsr.sort_values(by=['Group'])
        vals=dfarsr['Well ID'].values
        indexes=np.unique(vals, return_index=True)[1]
        myorder=[vals[index] for index in sorted(indexes)]

        plt.figure()
        sns.barplot(x='Well ID', y='awsr(min)', hue='Group', data=dfarsr, order=myorder)
        plt.savefig(to+'/'+'AWSRLines'+compound+label+'bar'+'StimPre'+'.tif', format='tif', bbox_inches='tight')
        plt.show()
        plt.figure()
        sns.catplot(x='Group', y='awsr(min)', data=dfarsr, kind='point')
        plt.savefig(to+'/'+'AWSRLines'+compound+label+'Point'+'StimPre'+'.tif', format='tif', bbox_inches='tight')
        plt.show()

        
        arsrled=ARSR(spikes_led, group_iN, durled, filt).values[0]
        dfarsr=pd.DataFrame(arsrled)
        dfarsr['awsr(min)']=dfarsr['awsr(min)'].apply(lambda x : np.log(x))

        #zscore well selection 

        dfarsr['z score']=np.abs(stats.zscore(dfarsr['awsr(min)']))

        #outlier well ids are defined in controls
        
       

        dfarsr=dfarsr[dfarsr['Well ID'].isin(SelectedWells)]

        dfarsr['Compound']=np.repeat(compound, len(dfarsr))
        dfarsr['Label']=np.repeat(label, len(dfarsr))
        dfarsr['Experiment']=np.repeat(ID, len(dfarsr))
        dfarsr['Light']=np.repeat('on', len(dfarsr))
        dfarsr['StimWells']=dfarsr['Well ID'].apply(lambda x: "on" if x is LED_IDs else 'off')
        df_ARSR=pd.concat([df_ARSR, dfarsr], axis=0).reset_index(drop=True)
        
        dfarsr=dfarsr.sort_values(by=['Group'])
        vals=dfarsr['Well ID'].values
        indexes=np.unique(vals, return_index=True)[1]
        myorder=[vals[index] for index in sorted(indexes)]

        plt.figure()
        sns.barplot(x='Well ID', y='awsr(min)', hue='Group', data=dfarsr, order=myorder)
        plt.savefig(to+'/'+'AWSRLines'+'bar'+'StimOn'+'.tif', format='tif', bbox_inches='tight')
        plt.show()
        plt.figure()
        sns.catplot(x='Group', y='awsr(min)', data=dfarsr, kind='point')
        plt.savefig(to+'/'+'AWSRLines'+'Point'+'StimOn'+'.tif', format='tif', bbox_inches='tight')
        plt.show()
        
        arsrpost=ARSR(spike_post, group_iN, durpost, filt).values[0]
        dfarsr=pd.DataFrame(arsrpost)
        dfarsr['awsr(min)']=dfarsr['awsr(min)'].apply(lambda x : np.log(x))

        #zscore well selection 

        dfarsr['z score']=np.abs(stats.zscore(dfarsr['awsr(min)']))

        #outlier well ids are defined in controls
        
        

        dfarsr=dfarsr[dfarsr['Well ID'].isin(SelectedWells)]

        dfarsr['Compound']=np.repeat(compound, len(dfarsr))
        dfarsr['Label']=np.repeat(label, len(dfarsr))
        dfarsr['Experiment']=np.repeat(ID, len(dfarsr))
        dfarsr['Light']=np.repeat('off_post', len(dfarsr))
        dfarsr['StimWells']=dfarsr['Well ID'].apply(lambda x: "on" if x in LED_IDs else 'off')
        df_ARSR=pd.concat([df_ARSR, dfarsr], axis=0).reset_index(drop=True)
        
        dfarsr=dfarsr.sort_values(by=['Group'])
        vals=dfarsr['Well ID'].values
        indexes=np.unique(vals, return_index=True)[1]
        myorder=[vals[index] for index in sorted(indexes)]

        plt.figure()
        sns.barplot(x='Well ID', y='awsr(min)', hue='Group', data=dfarsr, order=myorder)
        plt.savefig(to+'/'+'AWSRLines'+'bar'+'StimPost'+'.tif', format='tif', bbox_inches='tight')
        plt.show()
        plt.figure()
        sns.catplot(x='Group', y='awsr(min)', data=dfarsr, kind='point')
        plt.savefig(to+'/'+'AWSRLines'+'Point'+'StimPOst'+'.tif', format='tif', bbox_inches='tight')
        plt.show()


    else:
        
        seriesawsr=series[series['Dose Label']==label]
        #print(len(seriesawsr), 'lenforawsr')
        durbase=(end-start)/1000000
        
        arsrbase=ARSR(seriesawsr, group_iN, durbase, filt).values[0]
        dfarsr=pd.DataFrame(arsrbase)
        dfarsr['awsr(min)']=dfarsr['awsr(min)'].apply(lambda x : np.log10(x))

        #zscore well selection 

        dfarsr['z score']=np.abs(stats.zscore(dfarsr['awsr(min)']))

        #outlier well ids are defined in controls
        
        if len(SelectedWells)==0:

            SelectedWells=dfarsr[dfarsr['z score']<3]['Well ID'].values

        dfarsr=dfarsr[dfarsr['Well ID'].isin(SelectedWells)]

        dfarsr['Compound']=np.repeat(compound, len(dfarsr))
        dfarsr['Label']=np.repeat(label, len(dfarsr))
        dfarsr['Experiment']=np.repeat(ID, len(dfarsr))
        dfarsr['Light']=np.repeat('off', len(dfarsr))
        dfarsr['StimWells']=dfarsr['Well ID'].apply(lambda x: "off" if x in LED_IDs else 'off')
        
        
        df_ARSR=pd.concat([df_ARSR, dfarsr], axis=0).reset_index(drop=True)
        
        dfarsr=dfarsr.sort_values(by=['Group'])
        vals=dfarsr['Well ID'].values
        indexes=np.unique(vals, return_index=True)[1]
        myorder=[vals[index] for index in sorted(indexes)]
        
        if len(dfarsr)>0:

            plt.figure()
            sns.barplot(x='Well ID', y='awsr(min)', hue='Group', data=dfarsr, order=myorder)
            plt.savefig(to+'/'+'AWSRLines'+'bar'+'NOstim'+'.tif', format='tif', bbox_inches='tight')
            plt.show()
            plt.figure()
            sns.catplot(x='Group', y='awsr(min)', data=dfarsr, kind='point')
            plt.savefig(to+'/'+'AWSRLines'+'Point'+'NOstim'+'.tif', format='tif', bbox_inches='tight')
            plt.show()

        #print(df_ARSR[:10], 'Df_ARSR')

    return (df_ARSR, SelectedWells)
               


# In[ ]:





# In[ ]:





# In[112]:


def est_winsize(spiketrain):
    
    x=np.diff(spiketrain)
    x_max = max(x)
    x_min = min(x)
    N_MIN = 100   #Minimum number of bins (integer)
                #N_MIN must be more than 1 (N_MIN > 1).
    N_MAX = 5000  #Maximum number of bins (integer)
    N = range(N_MIN,N_MAX) # #of Bins
    N = np.array(N)
    D = (x_max-x_min)/N    #Bin size vector
    C = np.zeros(shape=(np.size(D),1))

    #Computation of the cost function
    for i in range(len(N)):
        edges = np.linspace(x_min,x_max,N[i]+1) # Bin edges
        ki = np.histogram(x,edges) # Count # of events in bins
        ki = ki[0]
        k = np.mean(ki) #Mean of event count
        v = sum((ki-k)**2)/N[i] #Variance of event count
        C[i] = (2*k-v)/((D[i])**2) #The cost Function
    #Optimal Bin Size Selection

    cmin = min(C)
    idx  = np.where(C==cmin)
    idx = int(idx[0])
    optD = D[idx]

    edges = np.linspace(x_min,x_max,N[idx]+1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(x,edges)
    
    
    fig = plt.figure()
    plt.plot(D,C,'.b',optD,cmin,'*r')
    plt.show()
    
    return optD
    


# In[113]:


def network_burst(df):
    
    df=df.sort_values(by='Start Timestamp [µs]')
    
    burst_times=df[' StartTimestamp [µs]'].values/1000
    
    diffs=np.diff(burst_times)
    
    nb_index=0
    
    
    for ind, diffs in enumerte(diffs):
        
        return ind
        
        
        
    


# In[ ]:





# In[ ]:





# In[114]:


def upload_join(paths, to):
    
    """for joined analysis of experiment it is better to upload data join and save in a new folder"""
    
    fnames=np.unique(np.array([p for  path in paths for p in  os.listdir(path)]))
    
    
    
    for fn in fnames:
        
        dat=pd.DataFrame()
        
        
        
        for path in paths:
            
            if fn in os.listdir(path):
                
                try:
                    updat=pd.read_excel(path+'\\'+fn, index_col=[0])
                    
                except:
                    
                    updat=pd.read_csv(path+'\\'+fn, index_col=[0])
                    
                dat=pd.concat([dat, updat], axis=0)
                
        dat.to_csv(to+'\\'+'joined'+fn.split('.')[0] +'.csv')
        

    


# In[115]:


###min max


# In[116]:


def MinMax(data, lim):
    """Lim is a featue range given by list [min, max]"""
    
    
    
    if len(np.unique(data))==1:
        
        
        data_scaled=data
    
    else:
        
        
        data_scaled=[(((lim[-1]-lim[0])*(x-np.min(data))/(np.max(data)-np.min(data)))+lim[0]) if np.unique(data).shape[0]!=1
                     else x
                     for x 
                     in data]
        
    
    
    return data_scaled
    
    


# In[117]:


def MinMaxall(dat, lim):
    """Lim is a featue range given by list [min, max]"""
    
    if dat.all()==0:
        data_scaled=dat
    
    else:
        data_scaled=[(((lim[-1]-lim[0])*(x-np.min(data))/(np.max(data)-np.min(data)))+lim[0]) if np.unique(data).shape[0]!=1
                     else x
                     for x 
                     in data for data in dat.T]
    
    return data_scaled
    
    


# In[118]:


def make_spider(df,columns, row, title, color):

    # number of variable
    categories=columns
    N = len(categories)

    # What will be the angle of each axis in the plot? (we divide the plot / number of variable)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    # Initialise the spider plot
    ax = plt.subplot(2, 2, row+1, polar=True)

    # If you want the first axis to be on top:
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    # Draw one axe per variable + add labels labels yet
    plt.xticks(angles[:-1], categories, color='grey', size=8)

    # Draw ylabels
    ax.set_rlabel_position(0)
    #plt.yticks([0,0.5,1], ["10","20","30"], color="grey", size=7)
    plt.ylim([-1.5,1.5])

    # Ind1
    values=df.loc[row].drop('Group').values.flatten().tolist()
    
    #print(values, 'val')
    values += values[:1]
    ax.plot(angles, values, color=color, linewidth=2, linestyle='solid')
    ax.fill(angles, values, color=color, alpha=0.4)

    # Add a title
    plt.title(title, size=11, color=color, y=1.1)


# In[119]:


import random


# In[121]:


class ComplexRadar():
    def __init__(self, fig, variables, ranges,
                 n_ordinate_levels=4):
        angles = np.arange(0, 360, 360./len(variables))

        axes = [fig.add_axes([0.1,0.1,0.9,0.9], polar=True,
                label = "axes{}".format(i))
    
                for i in range(len(variables))]
        
        
        l, text = axes[0].set_thetagrids(angles, 
                                         labels=[v[6:] if 'Change' in v else v for v in variables], fontsize=12)
        [txt.set_rotation(angle-90) for txt, angle 
             in zip(text, angles)]
        axes[0].tick_params(axis='both', which='major', pad=60)
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
        
        #print(data)
        sdata = _scale_data(data, self.ranges)
        
        ##print(self.angle, data, sdata[0])
        self.ax.plot(self.angle, np.r_[sdata, sdata[0]], label=g,  *args, **kw)
    def fill(self, data, *args, **kw):
        sdata = _scale_data(data, self.ranges)
        self.ax.fill(self.angle, np.r_[sdata, sdata[0]], *args, **kw)


# In[122]:


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


# In[ ]:





# In[ ]:





# In[ ]:





# In[165]:


from math import pi
from matplotlib import cm
from sklearn import preprocessing

from sklearn.preprocessing import MinMaxScaler


# In[356]:


def PSTH_Validation(df, df1, controlight, targetlight, to):
    
    exp=df['Experiment'].values[0]
    grp=df['Group'].values[0]
    dlabel=df['Label'].values[0]
    
   
    
    ##df is dataset containing all parameter features for each channel for each well
    ##PSTH s are for 50 ms windows for 10ms windows, first window is 40 window
    
    ##Validation holds if stimlation with targetlight and controllight light occured. 
    
    if [controlight, targetlight] or [targetlight, controlight] in np.unique(df['LightWave'].values):
        
        
       ##'Stim' column provides info on stimulation , on or off 
        
        ##plot control light vs df light for each cell group/experiment condition
        t1=df[(df['LightWave']==targetlight) & (df['Stim']=='On')]
        t2=df[(df['LightWave']==controlight) & (df['Stim']=='On')]
        t3=df1[(df['LightWave']==targetlight) & (df['Stim']=='On')]
        
        meltedt=pd.melt(t1, id_vars=['MFR(Hz)', 'Experiment', 'Compound', 'Stim', 'Label', 'Well ID',
       'Group', 'Rchange', 'Channel ID',  'PSTH20ratio', 'PSTH200ratio',
       'STTCpre', 'STTCpost', 'STTCStim', 'Lightcorr', 'LightLag', 'LightWave', 'Age'], var_name='PSTH', value_name='Peak')

        meltedt2=pd.melt(t2, id_vars=['MFR(Hz)', 'Experiment', 'Compound', 'Stim', 'Label', 'Well ID',
       'Group', 'Rchange', 'Channel ID',  'PSTH20ratio', 'PSTH200ratio',
       'STTCpre', 'STTCpost', 'STTCStim', 'Lightcorr', 'LightLag', 'LightWave', 'Age'], var_name='PSTH', value_name='Peak')
        
        melted3=pd.melt(t3, id_vars=['MFR(Hz)', 'Experiment', 'Compound', 'Stim', 'Label', 'Well ID',
       'Group', 'Rchange', 'Channel ID',
       'STTCpre', 'STTCpost', 'STTCStim', 'Lightcorr', 'LightLag', 'LightWave', 'Age', 'MFRpost'], var_name='PSTH', value_name='Peak')
        
        
        fig, ax=plt.subplots(figsize=(10, 8))
        sns.pointplot(data=meltedt, x='PSTH', y='Peak', join=False, color=targetlight, label=targetlight, ax=ax)
        sns.pointplot(data=meltedt2, x='PSTH', y='Peak', join=False, color=controlight, label=controlight, ax=ax)
        sns.pointplot(data=melted3, x='PSTH', y='Peak', join=False, color='c', label='early', ax=ax)
        
        psthcols=df.columns[df.columns.str.contains('PSTH')].values
        psthmax=max([float(txt.split('H')[1]) for txt in psthcols if txt.split('H')[1].isdigit()==True ]) #maxpsthwindow
        
        plt.title(grp, fontsize=20)



        plt.ylabel('Max Spike Count (50ms)', fontsize=20)

       

        ticklabels=np.arange(-50, psthmax, 50)
        #ticklabels[0]=ticklabels[0]+10

        plt.xticks(ticks=ax.get_xticks(), labels=ticklabels, fontsize=18)

        xspan=ax.get_xticks()

        plt.xlabel('PSTH', fontsize=20)

        


        plt.ylabel('Max Spike Count window', fontsize=20)
        plt.xlabel('Time(ms)')
        plt.yticks(fontsize=14)
        
        plt.ylim([0, 3])


        targetlight_patch =mpatches.Patch(color=targetlight, label='target')

        controllight_patch =mpatches.Patch(color=controlight, label='control')
        
        earlypatch=mpatches.Patch(color='c', label='First 5 trains')


        plt.legend(handles=[targetlight_patch, controllight_patch,  earlypatch])
        
        #plt.savefig(to+'/'+exp+grp+dlabel+'PSTHValsamescale'+'.png')
        
        plt.show()






        







# #to do for analysis
# 
# 1. Integrate PSTH early into  code 
# 2. Integrate PSTH_Plotting and PSTH_validation functions
# 3. Test 10ms windows for 1Hz stimulation 
# 4. Traslatable PSTH parameters for different drugs and cultures
#     1. NBQX treated peak is real peak for direct activation
#     2. Take integral for up to 200ms values and divide to Direct peak value
#     3. Or simply normalize with NBQX peak
#     

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


###plan is to get relative values normalized by nbqx value at this value is zero what to do 

###no, for the group, compute m, nbqx median and when take the relative values. 


# In[ ]:





# In[107]:


def normNBQX(m):
    
    
    joincols=m.columns[m.columns.str.contains('PSTH')==False]
    
    psthcols=m.columns[m.columns.str.contains('PSTH')==True]
    
    
    
    psths25=['PSTH'+str(float(i)) for i in np.arange(0, 35, 5)]
    
    normval=m[(m['Stim']=='On') & (m['Label']=='NBQX') ][psths25].mean().max()
    
    m[psthcols]=m[psthcols]/normval
    
    
    
    
    
    return m
    
    


# In[113]:


def func(data):
    plt.figure()
    
    sns.boxplot(data=data,  x='llightlag', y='lightlagvalue')
    plt.title(data['Group'].values[0])
    


# In[123]:


def Norm(df, delta):
    
    df['awsr(min)']=df['awsr(min)'].apply(lambda x:x-delta)
    
    return df

def Normalize(df):
    delta=df[df['Label']=='Vehicle']['awsr(min)'].values[0]-df[df['Label']=='Control']['awsr(min)'].values[0]
    dfnorm=df.groupby('Label').apply(lambda f: Norm(f, delta) if (f['Label'].any!='Control' or f['Label'].any!='Vehicle') else f).reset_index()
   
    return dfnorm
    


# In[124]:


def Ratio(f, Controlw, delta, dict_compares, df_from):
    
    s=f
    
    k=s.reset_index()['Label'].values[0]
    
    compare_ind=dict_compares[k]
    
    if compare_ind!='None':
    
        compare=df_from[df_from['Label']==compare_ind]['NEW_MFR'].values[0]
    
    
    
    
        f['awsr(min)']=f['awsr(min)'].apply(lambda x:((x-compare)/compare))
    
    return f

def Normalize(dfnorm, labels):
    
   
    
    doses=np.unique(dfnorm['Label'].values)
    
    
    compares=[]
    for dose in doses:
        
        ind=np.where(labels==dose)[0]
        
        if ind==0:
            compares.append(labels[ind].tolist()[0])
        else:
            compares.append(labels[ind-1].tolist()[0])
    dict_compares=dict(zip(doses, compares))
    
 
        
        
        
    df_from=dfnorm
    

    cont=dfnorm[dfnorm['Label']=='Control']['NEW_MFR'].values[0]
    
    delta=((dfnorm[dfnorm['Label']=='Vehicle']['NEW_MFR'].values[0]-cont)/cont) #'R of change for the Vehicle'
    
    dfnorm=dfnorm.groupby('Label').apply(lambda f: Ratio(f, cont, delta, dict_compares, df_from)).reset_index()
   
    return dfnorm
    


# In[126]:


def Ratio(f, Controlw, dict_compares, df_from):
    
    s=f
    
    k=s.reset_index()['Label'].values[0]
    
    compare_ind=dict_compares[k]
    
    if compare_ind!='None':
    
        compare=df_from[df_from['Label']==compare_ind]['NEW_MFR'].values[0]
    
    
    
    
        f['MFR(Hz)']=f['MFR(Hz)'].apply(lambda x:((x-compare)/compare))
    
    return f

def Normalize(dfnorm, labels):
    
   
    
    doses=np.unique(dfnorm['Label'].values)
    
    
    compares=[]
    for dose in doses:
        
        ind=np.where(labels==dose)[0]
        
        if ind==0:
            compares.append(labels[ind].tolist()[0])
        else:
            compares.append(labels[ind-1].tolist()[0])
    dict_compares=dict(zip(doses, compares))
    
 
        
        
        
    df_from=dfnorm
    

    cont=dfnorm[dfnorm['Label']=='Control']['NEW_MFR'].values[0]
    
    #delta=((dfnorm[dfnorm['Label']=='Vehicle']['NEW_MFR'].values[0]-cont)/cont) #'R of change for the Vehicle'
    
    dfnorm=dfnorm.groupby('Label').apply(lambda f: Ratio(f, cont, dict_compares, df_from)).reset_index()
   
    return dfnorm
    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


##Assignements function  test 


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




