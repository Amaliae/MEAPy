#!/usr/bin/env python
# coding: utf-8

# In[10]:


from statistics import mean
import json
import pickle
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from scipy.stats import wilcoxon
from scipy import stats

from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import LinearColorMapper, ColorBar

from datetime import date

import seaborn as sns
from scipy.stats import boxcox
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




import statsmodels.formula.api as smf

import random


# In[11]:


def find_folder(mainpath, folder, identifier=None, nonidentifier=None):
    """
    Search for a folder that contains 'csv' in its name in the specified directory.

    Parameters
    ----------
    directory : str
        The path of the directory to search in.

    Returns
    -------
    csv_folder : str or None
        The name of the folder containing 'csv' in its name, or None if no such folder is found.
    """
    # Iterate over the items in the specified directory
    
    return_path=''
    
    
    if folder in os.listdir(mainpath):
        
        item_path = os.path.join(mainpath, folder)
        
        if os.path.isdir(item_path):
            
            # Return the name of the folder if it contains 'csv' in its name
            
            if identifier==None or identifier in folder:

                return_path=item_path



            else:
                return_path=''

            if nonidentifier==None or nonidentifier not in folder:

                return_path=item_path

            else:
                return_path=''

            if len(return_path)>0:



                return return_path
        
            
    
    
    
    else:
        for item in os.listdir(mainpath):
        # Construct the full path of the item
            item_path = os.path.join(mainpath, item)
            

            # Check if the item is a directory and contains 'csv' in its name
            if folder in item:

                

                if os.path.isdir(item_path):
                # Return the name of the folder if it contains 'csv' in its name

                    if identifier==None or identifier in item:

                        return_path=item_path

                      

                    else:
                        return_path=''

                    if nonidentifier==None or nonidentifier not in item:

                        return_path=item_path

                    else:
                        return_path=''

                    if len(return_path)>0:



                        return return_path
            else:

                walk=list(os.walk(item_path))

                return_path=[os.path.join(path[0], string) for path in walk for string in path[1] if folder in string]

                if len(return_path)>0:



                    return return_path[0]



        
    # Return None if no folder with 'csv' in its name is found
    return None


# In[12]:


def Fold_Folders(mainpath):
    
    
    folder_path=mainpath+'\\'+'hdf5'
    
    metadatapath=mainpath+'\\'+'hdf5'
    
    save_to=mainpath+'\\'+'csv'
    
    
    if 'csv' not in os.listdir(mainpath):
        
        os.mkdir(save_to)
    
    
    
    dirs=os.listdir(mainpath)
    
    walk=list(os.walk(mainpath))
    
    walkpaths=[i for index in walk  for i in index if isinstance(i, str)==True]
    
    if 'hdf5' in dirs:
    
        folder_path=mainpath+'\\'+'hdf5'
        
    else:
        
        'Convert spikes and save in hdf5 folder'
        
    mwsdir=[i for i  in dirs if 'mws' in i ]
    
    
    if len(mwsdir)>0:
        
        metadatapath=[mainpath]
        
    elif len(walkpaths)>0:

        metadatapath=[]

        
        
        for wpth in walkpaths:
            
            dirswpth=os.listdir(wpth)
            
            mwsdirpth=[i for i  in dirswpth if 'mws' in i ]
            
            if len(mwsdirpth)>0:
                
                metadatapth=wpth
                metadatapath.append(metadatapth)
                

    

    
                
    else:
        print('no mwsin directory')
        
  
        
        
    return (folder_path, metadatapath, save_to)
 


# In[ ]:



        


# In[8]:


import neo

import quantities as pq

from elephant.spike_train_correlation import spike_time_tiling_coefficient


# In[9]:


def modify_cfg_tsfel(cfg_filemod):
    
    usespectral=['FFT mean coefficient', 'Fundamental frequency',
                 'Max power spectrum', 'Maximum frequency', 'Median frequency',
             'Power bandwidth', 'Spectral centroid', 'Spectral entropy', 
                 'Wavelet absolute mean', 'Wavelet energy', 'Wavelet entropy']
    usestatistical=list(cfg_file['statistical'].keys())
    usetemporal=['Area under the curve', 'Autocorrelation', 'Centroid', 
                 'Mean absolute diff', 'Mean diff', 'Median absolute diff',
                 'Median diff', 'Negative turning points', 'Neighbourhood peaks', 'Positive turning points', 
                 'Signal distance', 'Slope', 'Sum absolute diff', 'Zero crossing rate']
    
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
    
        if feat in usestatistical:

            cfg_filemod['temporal'][feat]['use']='yes'

        else:

             cfg_filemod['temporal'][feat]['use']='no'
                
    for feat in list(cfg_filemod['fractal'].keys()):
    
        cfg_filemod['fractal'][feat]['use']='no'
        
        
    return cfg_filemod
        
        


# In[10]:


def Channnel_Stat_Times(dfchannel):
    
    
    
    maxrate=max(dfchannel['spike rates'].values)
    
    minrate=min(dfchannel['spike rates'].values)
    
    medianrate=np.percentile(dfchannel['spike rates'].values, 90)
    
    maxtime=dfchannel[dfchannel['spike rates']==maxrate]['ttime'].values[-1]
    
    mintime=dfchannel[dfchannel['spike rates']==minrate]['ttime'].values[-1]
    
    mediantime=dfchannel[dfchannel['spike rates']==medianrate]['ttime'].values[-1]
    
   
    
    
    return pd.DataFrame(np.array([maxtime, mintime, mediantime]).reshape(-1, 3), columns=['MaxBin', 'MinBin', 'Per90Bin'])
    
    
    
    
    


# In[1]:


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
        


# In[12]:


def seminstateous_df(df, start, end, winsize, window):
    
    tmps=df['Timestamp [µs]'].values
    
    
    ar=np.zeros(int((end-0)/winsize)+2)
    
    #print(ar.shape, 'arshape')
    
   
    
    
    
    for tx, ty in zip(tmps[:-1], tmps[1:]):
        
        idxx=int(round(tx//winsize))
        
        idxy=int(round(ty//winsize))
        
       # print(idxx, idxy, tx, ty)
    
    
    
        delim=ty-tx
        
        rate=1000000/delim
        
        ar[idxx:idxy]=rate
        
        
    timeax=np.linspace(start, end+(2*winsize), len(ar[int(round(start/winsize)):]))-start
    
        
        
    dfroll=pd.DataFrame( ar[int(round(start/winsize)):])
    
    timeroll=pd.DataFrame(timeax)
    



    ###Rolling MEAN (other stat when mean)
    #for each channel individually
    
    
    arr=dfroll.rolling(window//winsize, closed='both', axis=0).mean().values[int(window//winsize)-1:, :]
    
    timearr=timeroll.rolling(window//winsize, closed='both', axis=0).mean().values[int(window//winsize)-1:, :]
    #print(arr.shape, timearr.shape)
    
    ms5dfrolled=pd.DataFrame({'1/ISI': arr.reshape(arr.shape[0]), 'Time': timearr.reshape(timearr.shape[0])})
    
    ms2df=pd.DataFrame({'1/ISI': ar[int(round(start/winsize)):], 'Time': timeax})
    
    
                       
        
    return ms5dfrolled, ms2df 


# In[ ]:





# In[13]:


def calc_semiinsta(df, LED, stamp, winsize, window):
    
    print(df['Dose Label'].unique())
    
    
    dlabel=df['Dose Label'].unique()[0]
    
    
    wellid=df['Well ID'].unique()[0]
    
    
    df_collect=pd.DataFrame()
    
    start, end=segment(df, stamp, dlabel, Duration=None)
    #ledarg=check_led(LED, stamp, df)
    
 
        
    lightlabels=led_label(LED, df, stamp)

    for index, name in enumerate(lightlabels):

        ledstart, ledend, leds, led_wells, light_label, light_dur=led_param(LED, df, stamp, [], index)

        if wellid in led_wells:


            for ld in leds:

                df5ms, df2ms=seminstateous_df(df, ld-50000, ld+1000000, winsize, window)

                df2ms['LighWaveform']=name
                df2ms['LED Start']=ld

                df_collect=pd.concat([df_collect, df2ms])

                
    return df_collect


# In[14]:


def Calc_tsfel_features(ch, tsfelfeatures, fs):
    
    feat=pd.DataFrame()
    
    chdat=ch['1/ISI'].values
    
    if len(chdat)>10:
    
        feat=tsfel.time_series_features_extractor(tsfelfeatures, chdat, fs=fs)
    
    
    return feat


# In[20]:




def return_tdsfel_featuresinuse():
    
    cfg_file = tsfel.get_features_by_domain()
    
    cfg_filemod=cfg_file.copy()

    usespectral=['FFT mean coefficient', 'Fundamental frequency',  'Max power spectrum', 'Maximum frequency', 'Median frequency',
                 'Power bandwidth', 'Spectral centroid', 'Spectral entropy', 'Wavelet absolute mean', 
                 'Wavelet energy', 'Wavelet entropy']
    usestatistical=list(cfg_file['statistical'].keys())
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

        if feat in usestatistical:

            cfg_filemod['temporal'][feat]['use']='yes'

        else:

             cfg_filemod['temporal'][feat]['use']='no' 
        
    for feat in list(cfg_filemod['fractal'].keys()):


        cfg_filemod['fractal'][feat]['use']='no'

    for feat in list(cfg_filemod['fractal'].keys()):

        cfg_filemod['fractal'][feat]['use']='no'

    return cfg_filemod


# In[15]:


def calc_semiinstanoled(df, LED, stamp, winsize, window, startmethod=''):
    
   
    dlabel=df['Dose Label'].unique()[0]
    
    if len(startmethod)>0:
        
        start=df[startmethod].values[0]-400000
        end=start+5000000
    else:
    

        start, end=segment(df, stamp, dlabel, Duration=None)
    
   
    df5ms, df2ms=seminstateous_df(df, start, end, winsize, window)


# In[16]:


def plot_raster_well(welldf, save_to):
    
    
    '''welldf is df with  merge col and spike labels and saves as html'''
    
    welllab=welldf['Well Label'].unique()[0]
    
    experiment=welldf['Experiment'].unique()[0]
    
    dlabel=welldf['Dose Label'].unique()[0]
    
    plt.figure()

    p1=figure(width=1000, height=1000, title=welllab+experiment+dlabel)
    output_file_path = save_to+'//'+welllab+experiment+dlabel+'Raster.html'
    output_file(output_file_path)
    

    chlabels=welldf['Channel Label'].unique()



    for indx, chlabel in enumerate(chlabels):

        spikesirreg=welldf[(welldf['Channel Label']==chlabel) & (welldf['_merge']==2)]['BurstSpikes'].values
        
        spikesburst=welldf[(welldf['Channel Label']==chlabel) & (welldf['_merge']==3)]['BurstSpikes'].values
        colors=return_bokeh_colormap('viridis')
        colors_range=np.linspace(0, 200, 256)
        div=np.diff(colors_range)[0]
        color_mapper = LinearColorMapper(palette = colors, low = 0, high=200)






        p1.rect(x=spikesirreg/1000,  y=np.ones(len(spikesirreg))*chlabel, width=2, height=0.5,
                fill_color='blue', line_color="blue")
        
        p1.rect(x=spikesburst/1000,  y=np.ones(len(spikesburst))*chlabel, width=2, height=0.5,
                fill_color='red', line_color="red")
       

    show(p1)
    save(p1, filename=output_file_path)
    


# In[17]:


def MFRSIrreg(diffwhereirreg, spkwhereirreg):
    
    #alliregreg=[]

    currireg=[]

    mfrs=[]
    
    
    


    for i, j in zip(diffwhereirreg, spkwhereirreg[:-1]):

        if i==1:

            currireg.append(j)

        else:
            
            #alliregreg.append(1000000/np.diff(np.array(currireg)))
            
            if len(currireg)>1:

                mfrs.append(len(currireg)/((currireg[-1]-currireg[0])/1000000))

            currireg=[]
            pass
        
    
        
    return mfrs#, alliregreg
       


# In[18]:


def Calculate_Irregular_Rates(mergeddfbychannel):
    
    '''in merge 3 is both, 2 is left only, so non bursting spike'''
    
    
    
    mergeddfbychannel=mergeddfbychannel.sort_values(by='BurstSpikes')

    spk=mergeddfbychannel['BurstSpikes'].values
    
    

    spkisi=1000000/np.diff(spk)
    
    labelsarr=np.ones(len(spkisi))

    spktype=mergeddfbychannel['_merge'].values

    whereirreg=np.where(spktype==2)[0]

    spkwhereirreg=spk[whereirreg]

    diffwhereirreg=np.diff(whereirreg)

    whereirregcon=whereirreg[np.where(np.diff(whereirreg)==1)[0]]


    spkirregisi=spkisi[whereirregcon]
    
    labelsarr[whereirregcon]=0

    spkmfrs=MFRSIrreg(diffwhereirreg, spkwhereirreg)
    
    new_df=pd.DataFrame({'1/ISI':spkisi, 'Irregular':labelsarr})
    
    return new_df    


# In[2]:


def load_object(filepath):
    
    with open(filepath, 'rb') as inp:
        
        load_obj= pickle.load(inp)

    return load_obj


# In[19]:


def modify_burstdata(dat):
    
    df=pd.DataFrame({'BurstSpikes': dat['Spikes'].values[0]})
    
    return df


# In[20]:


def modify_networkburstdata(dat):
    
    ar=dat['Spikes per channel'].values[0]
    chl=dat['Channel Labels'].values[0]
    
    nbt=[i[0] for i in ar]
    
    
    df=pd.DataFrame({'Timestamp [µs]': nbt, 'Channel Label':chl})
    
    return df


# In[21]:


def led_seg_old(led, df, stamp):
    
    file=df['File'].unique()[0]
    
    pos=file.find('h5')+2
    
    
    
    key1=file[:pos]
    key2=df['Dose Label'].unique()[0]
    
    dose_start, dose_end=segment(df, stamp, key2)
    
    times=led[key1]['Time'][0]
    
    if len(times[(times>dose_start) & (times<dose_end)])>0:
    
    
    
        start=times[(times>dose_start) & (times<dose_end)][0]
        end=times[(times>dose_start) & (times<dose_end)][-1]
    
    
    print(start, end)
    
    
    return (start, end)
    


# In[ ]:





# In[ ]:





# http://fmwww.bc.edu/repec/bocode/t/transint.html text on transformations (wo author info) that makes sense to me

# In[22]:


def find_peak_slope(xx, yy):
    
    
    peak=max(yy)
    
    
    
    

    peakindex=np.argmax(yy)
    
    riseslope=peak/peakindex
    
    fallslope=-1*peak/peak
    
    

    peak50=np.where(yy<=peak*0.5)[0]
    
   
    
    if len(peak50)>0:
        
      

        peak50indexes=peak50[peak50<peakindex]
        
        peak50findexes=peak50[peak50>peakindex]
        
        if len(peak50indexes)>1:
           
            
            peak50index=peak50indexes[-1]
            
            riseslope=(yy[peakindex]-yy[peak50index])/(peakindex- peak50index)
            
        if len(peak50findexes)>1:
            
         
            
            peak50findex=peak50findexes[0]
            
            fallslope=(yy[peak50findex]-yy[peakindex])/(peak50findex-peakindex)
        
        
   #plt.plot(xx, (peak*0.5)+xx*riseslope, label='rise')
   # plt.plot(xx, yy, label='curve')
    
    #plt.plot(xx, peak+xx*fallslope, label='fall')''
    
    
    
    return (riseslope, fallslope)
    


# In[13]:


def burstpolyfit(x, y, degree):
    #degree is 8 
    args={}

    coeff= np.polyfit(x, y, degree)
    
    coeff2=np.polyfit(x, y, 1)
    
    
    #print(coeff, coeff)
   
    
    func=np.poly1d(coeff)
    func2=np.poly1d(coeff2)
    ynew=func(x)
    ynew1=func2(x)
    
    
    
   # plt.figure()


    #plt.plot(x, y, '.')
    #plt.plot(x, ynew1, '--', label='Polyfit')
    
    #plt.plot(x, ynew, '-')
    
    args['coeff']=coeff[::-1]
    
    
    args['x']=x
    
    args['y']=y
    
    args['yhat']=ynew
    
    
    return args

def zoom_burst(b, params):
    
    firstN=params[0]
    winsize=params[1]
    degree=params[2]
    instat=params[3]
    windowsize=params[4]
    
    args=fit_curve(b, firstN, winsize, degree, insta=instat)
    
    slopes=find_peak_slope(args['x'], args['yhat'])
    
    area=area_under_curve(args['x'], args['yhat'], windowsize)
    
    args['Rise']=slopes[0]
    
    args['Fall']=slopes[1]
    
    args['Area']=area
    
    return args
    
    
    
def area_under_curve(xx, yy, Nspikes):
    
    xind=(xx>0)& (xx<Nspikes)
    
    area=np.trapz(yy[xind], xx[xind])
    
    #plt.plot(xx, yy)
    
    
    #plt.fill_between(xx[xind], y1=yy[xind])
    
    
    
    #plt.legend()
    
    #plt.show()
    
    return area


def fit_curve(b, firstN, winsize, degree, insta=True):
    
    mb=modify_b(b)
    
    prey=mb[:firstN]
    
    if insta==True:
        
        y=1000000/np.diff(prey)
        
        x=np.arange(0, len(y), 1)
    
        args=burstpolyfit(x, y, degree)
        
    else:
        
        y=bin(prey, prey[0], prey[-1], winsize)
        
        x=np.linspace(prey[0],  prey[0]+(len(y)-1)*winsize, len(y))/1000000
        
        args=burstpolyfit(x, y, degree)
        
        
    
    return args
    
    


# In[25]:





# In[26]:





# In[27]:


def burst_rate(df, dictdur):
    
    expid=df['Experiment'].unique()[0].find('iN')
    exp=df['Experiment'].unique()[0][expid:expid+4]
    
    dur=dictdur[exp]
    
    burstr=len(df['Timestamp [µs]'].values)/dur
    
    
    
    return burstr


# In[28]:


def relative_c(df, tosub, hue_order, Params, groubby, to, save):
    
    
    
    print(df[:3], 'dgdg')
    
    compound=df['Compound ID'].unique()[0]
    
    print(compound, 'compound')
    
   
   
    
    if compound!='Control' and compound!='No Group':
        
    

        ide=df['Experiment']
        
        
        print(ide, 'ideee???')
        
       
        print(hue_order, 'hue_order')

       

        labels=np.unique(df['Dose Label'].values)
        
        print(labels, 'labels')

        newcols=['Change'+col for col in Params]
        dfnew=pd.DataFrame(columns=newcols)

        for label in labels:

            if label!='Control':

                ind=np.where(np.array(hue_order)==label)[0]
                
                

                print(ind, label, 'inlopdjgh')
                
                print(np.array(hue_order)[int(ind-1)], 'ind-1')

                dfpre=df[df['Dose Label']==hue_order[int(ind-1)]].sort_values(by='Channel ID')

                dfpost=df[df['Dose Label']==hue_order[int(ind)]].sort_values(by='Channel ID')

                x=dfpre[Params].values
                y=dfpost[Params].values

                r=rchange(x, y)
                
                print(len(r), 'length', len(dfpre), len(dfpost), label)



                dfpost[newcols]=r

               

                dfnew=pd.concat([dfnew, dfpost], axis=0)
                
        if save=='True':
                
        

            dfnew.to_csv(to+'/Files'+'/Change of'+Param+'.csv') 
            
            
        return dfnew## changed this ####


# In[29]:


def instrates(df, start, stop, winsize, param):
  
    
    spikes=df['Timestamp [µs]'].values
    
    
    chrates=bin(spikes, start, stop, winsize)
    
    time=np.arange(1, len(chrates)+1, 1)
    
    return pd.DataFrame({param+' '+'rates':np.array(chrates), 'btime':time})
    
   


# In[30]:


def instratesTrueTime(df, start, stop, winsize, param):
  
    
    spikes=df['Timestamp [µs]'].values
    
    
    chrates=bin(spikes, start, stop, winsize)
    
    truetime=np.arange(1, len(chrates)+1, 1)*winsize+start
    
    time=np.arange(1, len(chrates)+1, 1)
    
    return pd.DataFrame({param+' '+'rates':np.array(chrates), 'btime':time, 'ttime':truetime})
    
   


# In[31]:


def relative_here2(dburst, params, tosub, hue_order):
    
    
    
    ss=dburst.groupby([ 'Compound ID', 'Group'])['Dose Label'].value_counts().min()
    ndburstdir=dburst.groupby(['Compound ID', 'Group', 'Dose Label']).apply(lambda x: 
                                                x.iloc[random.sample(range(0, len(x)), ss), :]).reset_index(drop=True)
    #params=['Duration', 'Spike Count', 'Max ISI', 'Min ISI', 'Mean ISI',
        #'Burstiness', 'FT Spikes',
        #'Surprise']
    
    
    ndburst=ndburstdir.groupby(['Compound ID', 'Group']).apply(lambda x:
                                                                        relative_c(x, tosub, hue_order,
                                                                    params, 'Group',  tosub, 'False').reset_index(drop=True))
        
    
 
        
        
    return (ndburst.reset_index(drop=True), ndburstdir)


# In[32]:


def percentile_rates(rates, param, ps):
    
    """transform is a dict with the name of transform is ilst  with a first val and sec val param if needed"""
    
    df=pd.DataFrame()
    
    vals=rates[param].values
    
    for index, percent in enumerate(ps):
        
        
        if index!=0:
            
            rangelow=np.percentile(vals, ps[index-1])
            
            rangehigh=np.percentile(vals, percent)
            
            rangeval=vals[(vals>=rangelow) & (vals<rangehigh)]
            
          
            df[str(ps[index-1])+str(percent)]=rangeval
                
                
                
            
    return df
            
        
        
    
    
    
    
    
    


# In[33]:


def rate_prop(dfall, percsval, param):
    
    """dfall is a grouped by object, percs is a list of percentiles of which proportions are caclulated"""
    
    
    rts=dfall[param].values
    
    
    df=pd.DataFrame()
    
   
    
    
    
    for pi in range(1, len(percsval)):
        
        rangeval=len(rts[(rts<percsval[pi]) & (rts>=percsval[pi-1])])
        
        
        
        
        df[str(int(percsval[pi-1]))+str(int(percsval[pi]))]=[rangeval/len(rts)]
        
    rangeval=len(rts[rts<percsval[0]])
    
    df[str(int(percsval[0]))+str(int(percsval[0]))]=[rangeval/len(rts)]
    
    rangeval=len(rts[rts>percsval[-1]])
    
    df[str(int(percsval[-1]))+str(int(percsval[-1]))]=[rangeval/len(rts)]
        
    
        
        
        
        
    return df
            


# In[34]:


def rate_prop0(dfall, percs, param):
    
    """dfall is a grouped by object, percs is a list of percentiles of which proportions are caclulated"""
    
    
    rts=dfall[param].values
    
    
    df=pd.DataFrame()
    
    percsval=[np.percentile(rts, pc) for pc in percs]
    
    for pi in range(1, len(percsval)):
        
        rangeval=len(rts[(rts<percsval[pi]) & (rts>percsval[pi-1])])
        
        
        df['Prop'+str(percsval[pi-1])+str(percsval[pi])]=rangeval/len(rts)
        
        
        
        
    return df
            


# In[35]:


def transform_param(rangeval, transform):
    
    if transform[0]=='boxcox':
        
        trasformedvar=pd.Series(boxcox(rangeval, lmbda=transform[0]))
                
    elif transform[0]=='reciprocal root':
        
        trasformedvar=-1*rangeval**-0.5
        
    else:
        
        transformedvar=rangevar
        
    return rangevar
        
        
        
        


# In[36]:


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


# In[37]:


def load_exportedext(path, identifier, identifier2, identifier3, newpharm, Comps):
    
    
    
    dirs = os.listdir(path)
    
    df_bursts=0
    df_data=0
    
    
    filename=[]
    
    i=0
    j=0
    
   
    
    for f in dirs:
        
        #print(f)
        #print(f[:-10])
       
        
        if re.search(identifier, f) and re.search(identifier2, f) and re.search(identifier3, f):
            
            if re.search('.csv$', f):
            
               
            
                read_func=pd.read_csv
            
            elif re.search('.pkl$', f):
            
               

                read_func=pd.read_pickle
            
            else:
                read_func=pd.read_xlsx
            
            
        
        
            
            if re.search(identifier, f):
                
                
                
                if i==0:
                    
                    df_data=read_func(path+'/'+f)
                    
                    if 'Unnamed: 0' in df_data.columns:
                        
                        df_data.drop(['Unnamed: 0'], axis=1, inplace=True)
                    
                    df_data['File']=f
                    #df_data['Experiment']=df_data['Experiment'].apply(lambda x : f[:-10])
                    #df_data['Dose Label']=df_data['Dose Label'].apply(lambda x : str(x))

                    i=i+1

                else:
                    to_c=read_func(path+'/'+f)
                    to_c['File']=f
                    
                    if 'Unnamed: 0' in to_c.columns:
                        to_c.drop(['Unnamed: 0'], axis=1, inplace=True)
                    #to_c['Experiment']=to_c['Experiment'].apply(lambda x : f[:-10])
                    #to_c['Dose Label']=to_c['Dose Label'].apply(lambda x: str(x))
                    df_data=pd.concat([df_data, to_c])
                    
                
            filename.append(f)
            
          
    
    try:
        A=dict(zip(np.unique(df_data['Channel Label'].values), np.arange(0, 12, 1)))
        df_data['True Label']=df_data['Channel Label']
        df_data['Channel Label']=df_data['Channel Label'].apply(lambda x : A[x])
        
    except: 
        print('No spike files')
 
    try:
        
        df_data['Dose Label']=df_data['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)
        
        
        if newpharm==True:
            
            
            df_data['Dose Label']=df_data['Compound ID']
            
            
            df_data['Compound ID']=df_data['Well ID'].apply(lambda x : Comps[int(x)])
            
            #print(np.unique(df_data['Compound ID'].values), 'inloadexpoetedcompaounds')
            
            #print(np.unique(df_data['Compound ID'].values), 'after')
            
    except: 
        print('No spike files')
    
        #df_bursts['Dose Label']=df_bursts['Dose Label'].apply(lambda x : 'Control' if x=='0' else x)
        
   

    columnsburst=df_data.columns
        
        #print(columnsburst, 'cols')

    columnsburst=['Timestamp [µs]' if column=='Start timestamp [µs]' else column for column in columnsburst]

        #print(columnsburst, 'columnsb')

    df_data.columns=columnsburst
        
    return (df_data, filename)


# In[39]:


def load_data_fromlist(dpaths, names, Type, newpharm, Comps, identifier1, identifier2):
       
       
   dfall=pd.DataFrame()
   
   
   for path, dname in zip(dpaths, names):
       
       if Type=='Spike':
           
           try:
               f=load_exported(path, identifier1, identifier2, newpharm, Comps)[0]
               
               

          
               f.name=dname

               dfall=pd.concat([dfall, f], axis=0)

           except:
               
               print('No spike data')
               
       if Type=='Burst':
           
           try:
               f=load_exported(path, identifier1, identifier2, newpharm, Comps)[1]

          
               f.name=dname

               dfall=pd.concat([dfall, f], axis=0)

           except:
               
               print('No burst data')

   return dfall


# In[40]:


def load_meta_fromlist(dpaths, names):
    
    Groupall={}
    stampall={}
    LEDall={}
    Bookall={}
    
    
    for dpath, name in zip (dpaths, names):
        
        stamp, LED, Book, Group_iN=load_meta(dpath)
        
        Groupall[name]=Group_iN
        Bookall[name]=Book
        
        stampall.update(stamp)
        LEDall.update(LED)
        
        
    return stampall, LEDall, Bookall, Groupall
            


# In[ ]:





# In[41]:


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


# In[42]:


def get_best_distribution(data):
    
    ###TODO copied from somwhere , check the code!!
    
   
    
    import scipy.stats as st
    
    dist_names = ["norm", "beta",  "exponweib", "weibull_max", "weibull_min", "pareto", 
                  "genextreme", "rayleigh", "lognorm"]
    dist_results = []
    params = {}
    for dist_name in dist_names:
        dist = getattr(st, dist_name)
        param = dist.fit(data)

        params[dist_name] = param
        # Applying the Kolmogorov-Smirnov test
        D, p = st.kstest(data, dist_name, args=param)
        print("p value for "+dist_name+" = "+str(p))
        dist_results.append((dist_name, p))

    # select the best fitted distribution
    best_dist, best_p = (max(dist_results, key=lambda item: item[1]))
    # store the name of the best fit and its p value

    print("Best fitting distribution: "+str(best_dist))
    print("Best p value: "+ str(best_p))
    print("Parameters for the best fit: "+ str(params[best_dist]))

    return (best_dist, best_p, params[best_dist])


# In[43]:


def load_object(filepath):
    
    with open(filepath, 'rb') as inp:
        
        load_obj= pickle.load(inp)

    return load_obj


# In[44]:


def param_DistStatDPL(df, param):
    
    """take values by groupby and calculation  """
    
    df['Sum'+' '+param]=sum(df[param].values)
    
    df['Q01'+' '+param]=np.quantile(df[param].values, 0.1)
    df['Q1'+' '+param]=np.quantile(df[param].values, 0.25)
    df['Median'+' '+param]=np.quantile(df[param].values, 0.5)
    df['Q3'+' '+param]=np.quantile(df[param].values, 0.75)
    df['Q4'+' '+param]=np.quantile(df[param].values, 0.95)
    df['Max'+' '+param]=np.quantile(df[param].values, 1)
    df['Mean'+' '+param]=np.mean(df[param].values)
    df['STD'+' '+param]=np.std(df[param].values)
    df['CV'+' '+param]=stats.variation(df[param].values, ddof=1, nan_policy='omit')
    
    return df.iloc[0, :]


# In[45]:


def param_DistStat(df, param):
    
    """take values by groupby and calculation, modified not to have a duplicate """
    dfnew={}
    
    dfnew['Sum'+' '+param]=sum(df[param].values)
    
    dfnew['Q01'+' '+param]=np.quantile(df[param].values, 0.1)
    dfnew['Q1'+' '+param]=np.quantile(df[param].values, 0.25)
    dfnew['Median'+' '+param]=np.quantile(df[param].values, 0.5)
    dfnew['Q3'+' '+param]=np.quantile(df[param].values, 0.75)
    dfnew['Q4'+' '+param]=np.quantile(df[param].values, 0.95)
    dfnew['Max'+' '+param]=np.quantile(df[param].values, 1)
    dfnew['Mean'+' '+param]=np.mean(df[param].values)
    dfnew['STD'+' '+param]=np.std(df[param].values)
    dfnew['CV'+' '+param]=stats.variation(df[param].values, ddof=1, nan_policy='omit')
    
    return pd.DataFrame(dfnew, index=[0])


# In[46]:


def plot_config(plotparamstochange, values, nrows, ncols):
    
    
    plotparam={}
    plotparam['style']='ggplot'
    plotparam['context']='paper'
    
    plotparam['sharey']=True
    
    plotparam['figsize']=(16, 6)

    plotparam['dpi']=600

    plotparam['ticklabelsize']=14
    plotparam['labelsize']=16
    plotparam['formats']=['.svg']
    plotparam['colours']=[]

    plotparam['axnrows']=2
    plotparam['type']=['mfr', 'changemfr']

    plotparam["font.family"] = "Arial"

    for param, value in zip(plotparamstochange, values):
        
        print(plotparam[param], value)
        
        plotparam[param]=value
        
        
    
    sns.set_context(plotparam['context'])

    plt.rc('xtick', labelsize=plotparam['ticklabelsize']) 
    plt.rc('ytick', labelsize=plotparam['ticklabelsize'])
    plt.rc('axes', labelsize=plotparam['labelsize'])
    plt.rcParams["font.family"] =plotparam["font.family"]# fontsize of the x and y labels
    
   

            



    fig, axs = plt.subplots(nrows, ncols, figsize=plotparam['figsize'],  sharey=plotparam['sharey'])
    #constrained_layout=True
    
    
    return fig, axs, plotparam


# In[47]:


def plot_how(plotparamstochange, value, plotype, dts, compare_order, hue_order, colour_dict, AgeParam, nrows, annotit, to):
    
    """style your plots dependet of what type of plot
    
    return plotarg"""
    
    stat_dict={'Labels':{}, 'Groups':{}, 'Age':{}}
    
    
    if plotype=='Age':
        
        data=dts[0]
        
        y=dts[2]
        
        ncols=len(compare_order)
        

        ageorder=data.sort_values(by=AgeParam)[AgeParam].unique().astype('str')

        data[AgeParam]=data[AgeParam].astype('str')


        
        ticksize=14
        labelsize=16


        scale = 1 # scaling factor for the plot

        fig_width = (ncols *scale*4)
        fig_height = 4

        fig, axs, pp=plot_config(['figsize'], [(fig_width, fig_height)],  1, 1)

        #ageorder=dat.sort_values(by='Age')['Age'].unique().astype('str')

        #dat['Age']=data['Age'].apply(lambda x: str(x))
        testres=stat_testpub(data, AgeParam, 'Group', 'indep')


        signpvalues=data.groupby(AgeParam).apply(lambda x: statistical_test(x, 'Group', y, 'Ch', stats.mannwhitneyu, False,  [])[1])
        signpairs=data.groupby(AgeParam).apply(lambda x: statistical_test(x, 'Group', y, 'Ch', stats.mannwhitneyu, False,  [])[2])
        dit=list(signpairs.to_dict().items())

        pvaldict=list(signpvalues.to_dict().items())
        pairs=[((pair[0],  p[0]), (pair[0], p[1]))  for pair in dit for p in pair[1]]
        pvalues=[pvaldict[i][1][p] for i in range(len(dit)) for p in dit[i][1]]
        
        stat_dict[AgeParam]=dict(zip(pairs, pvalues))   
        #pairs=[pair for pair in pairs for l in list(signpairs[signpairs.index==int(pair[0][0])].values[0].keys()) if pair[0][1] and pair[1][1] in l]

        sns.boxplot(data=data[(data['Label']=='Control')], x=AgeParam,  y=y, ax=axs, hue='Group',  width=0.5, order=ageorder,
                               hue_order=compare_order, palette=colour_dict.values(), boxprops=dict(alpha=.7), showfliers=False)

        annot = Annotator(axs, pairs, data=data[(data['Label']=='Control')],  x=AgeParam,  y=y, hue='Group', order=ageorder, 
                                hue_order=compare_order, palette=colour_dict.values(),  boxprops=dict(alpha=.7), showfliers=False)#, order=order, hue_order=hue_order)

        #annot.configure(test=test1, text_format='star', loc='outside', comparisons_correction="bonferroni")
        #annot.apply_test()

        annot.configure(loc='outside')


        annot.set_pvalues(pvalues)


        annot.annotate()

        #annot.configure(loc='outside')

        axs.set_xlabel(AgeParam, fontsize=labelsize)

                #b.axes.set_title("Title",fontsize=50)

        axs.tick_params(labelsize=ticksize)
        
        fig.legend(bbox_to_anchor=(1.07, 0.9), fontsize=labelsize)
        
        axs.get_legend().remove()
            
        
        plt.savefig(to+'/'+AgeParam+y+'.svg', bbox_inches='tight')






        plt.show()   

    
    
    if plotype=='MFRs':
        
        
        ncols=len(compare_order)
        
        data=dts[0]
        x=dts[1]
        y=dts[2]
        ticksize=10
        labelsize=12
       
        
        scale = 1 # scaling factor for the plot
        
        fig_width = (ncols *scale*2.5)
        fig_height = 2.5

        fig, axs, pp=plot_config(['figsize'], [(fig_width, fig_height)],  nrows, ncols)

        for ax, g in zip(axs[:], compare_order):
            
            color=colour_dict[g]

            sns.boxplot(data=data[data['Group']==g], x='Label',  y=dts[2], ax=ax, color=color, width=0.5, order=hue_order)
            
            
            ax.set_xlabel(g, fontsize=labelsize)
            
            #b.axes.set_title("Title",fontsize=50)
           
            ax.tick_params(labelsize=ticksize)
            
            
        print(x, y)
            
        plt.savefig(to+'/'+x+y+'Label'+'.svg', bbox_inches='tight')

      

        

   
        
        ncols=len(hue_order)
        
       
        scale = 1# scaling factor for the plot
       
        
        
        
       
        
        handles=[mpatches.Patch(color=v, label=k) for k, v in colour_dict.items()]
        
       
        
        fig_width = (ncols*scale*2.5)
        fig_height = 2.5
        
       

        fig, axs, pp=plot_config(['figsize', 'sharey'], [(fig_width, fig_height), False], nrows, ncols)
        
        axisindex=-1
      

       

        for ax, g in zip(axs[:], hue_order):
            
            axisindex=axisindex+1
            
          
            sns.boxplot(data=data[data['Label']==g], x='Group',  y=dts[2], ax=ax, width=0.5,  
                        dodge=False, palette=colour_dict.values(), order=compare_order)
            
            #if g=='Control':
                
                
                
            pair, test1=stat_testpub(data[data['Label']==g], 'Group', 'Label', 'indep')

            annot = Annotator(ax, pair, data=data[data['Label']==g], x='Group',  y=dts[2],
                       order=compare_order, palette=colour_dict.values(),  boxprops=dict(alpha=.7))#, order=order, hue_order=hue_order)


            annot.configure(test=test1, text_format='star', loc='outside')
            annot.apply_test()

            if annotit==True:
                _, test_results =annot.annotate()



                stat_results = [(result.data.stat_value, result.data.pvalue) for result in test_results]

                stat_groups=[(result.__dict__['structs'][0]['group'], result.__dict__['structs'][1]['group'])
                       for result in test_results]


                stat_dict['Groups'][g]=dict(zip(stat_groups, stat_results))    

                
                
                
                
                
                
                
                
            ax.set_xlabel(g, fontsize=labelsize)
            ax.set_xticks([])
            #b.axes.set_title("Title",fontsize=50)
           
            ax.tick_params(labelsize=ticksize)
            
            if axisindex>0:
                ax.set_ylabel(' ')
            
            
            #ax.get_legend().remove()
            
         
            
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.9, 1), fontsize=labelsize)
        
        print(x, y, handles)
            
        plt.savefig(to+'/'+x+y+'.svg', bbox_inches='tight')

       

           

            
    if plotype=='MFRchange':
        
        ticksize=10
        labelsize=12
        
        data=dts[0]
        
        x=dts[1]
        y=dts[2]
        
        ncols=len(compare_order)
        
        nrows=1
        
       
        scale = 1# scaling factor for the plot
       
        
        fig_width = (ncols*scale*2.5)
        fig_height = 2.5
        
        print(fig_width)

        fig, axs, pp=plot_config(['figsize', 'sharey'], [(fig_width, fig_height), False], nrows, ncols)

        axindex=-1
        
        for ax, g in zip(axs[:], compare_order):
            
            axindex=axindex+1
            
            color=colour_dict[g]
            
            hue_orderplot=np.delete(np.array(hue_order), int(np.where(np.array(hue_order)=='Control')[0]))
            
            
            pair, test1=stat_testpub(data[data['Group']==g], 'Label', 'Group',  'dep')
            
            #print(pair)

            sns.boxplot(data=data[data['Group']==g], x='Label', y=dts[2], ax=ax, width=0.5, 
                        showfliers=not(True==annotit), order=hue_orderplot, color=color)
            
            pair=[p for p in pair
                  if (np.where(np.array(hue_orderplot)==p[0])[0]-np.where(np.array(hue_orderplot)==p[1])[0])==-1]
            
            
            annot = Annotator(ax, pair, data=data[data['Group']==g], x='Label', y=dts[2],
                           order=hue_orderplot)#, order=order, hue_order=hue_order)


            annot.configure(test=test1, text_format='star', loc='outside')
            annot.apply_test()
            
            if annotit==True:
                _, test_results =annot.annotate()
                
                
                
                stat_results = [(result.data.stat_value, result.data.pvalue) for result in test_results]
                
                stat_groups=[(result.__dict__['structs'][0]['group'], result.__dict__['structs'][1]['group'])
                       for result in test_results]
                
               
                stat_dict['Labels'][g]=dict(zip(stat_groups, stat_results))    
                
            else:
                
                ax.set_ylim([-2, 10])

            ax.set_xlabel(g, fontsize=labelsize)
            
            

            #b.axes.set_title("Title",fontsize=50)
           
            ax.tick_params(labelsize=ticksize)
            
            ylabel=gen_y_label(dts[2])
            
            ax.set_ylabel(ylabel, fontsize=labelsize)
            
            if axindex>0:
                 ax.set_ylabel(' ')
            
            
           
            
           
            
            #ax.get_legend().remove()
            
        #plt.legend(loc='upper left', bbox_to_anchor=(1.03, 1), fontsize=labelsize)
        
        plt.savefig(to+'/'+x+y+str(annotit)+'.svg', bbox_inches='tight')
       
        hue_orderplot=np.delete(np.array(hue_order), int(np.where(np.array(hue_order)=='Control')[0]))
    
        
        
        ncols=len(hue_orderplot)
        
        nrows=1
        
       
        scale = 1# scaling factor for the plot
       
        
        fig_width = (ncols*scale*3)
        fig_height = 3
        
        print(fig_width)
        
        
        handles=[mpatches.Patch(color=v, label=k) for k, v in colour_dict.items()]
      
        fig, axs, pp=plot_config(['figsize', 'sharey'], [(fig_width, fig_height), False], nrows, ncols)

        print(handles, 'handles')

        axisindex=-1  
        
        for ax, g in zip(axs[:], hue_orderplot):
            
            axisindex=axisindex+1
            
            
            
            pair, test1=stat_testpub(data[data['Label']==g], 'Group', 'Label', 'indep')
            
            
            
            #print(pair)

            sns.boxplot(data=data[data['Label']==g], x='Group',  y=dts[2], ax=ax, width=0.5, 
                        showfliers=not(True==annotit), order=compare_order, palette=colour_dict.values(), boxprops=dict(alpha=.7))
            
            #pair=[p for p in pair
                  #if (np.where(np.array(hue_orderplot)==p[0][1])[0]-np.where(np.array(hue_orderplot)==p[1][1])[0])==-1]
            
            #print(pair)
            annot = Annotator(ax, pair, data=data[data['Label']==g], x='Group',  y=dts[2],
                           order=compare_order, palette=colour_dict.values(),  boxprops=dict(alpha=.7))#, order=order, hue_order=hue_order)


            annot.configure(test=test1, text_format='star', loc='outside', comparisons_correction="bonferroni")
            
            annot.apply_test()
           
            
            if annotit==True:
                
                _, test_results =annot.annotate()
                
                
                
                stat_results = [(result.data.stat_value, result.data.pvalue) for result in test_results]
                
                stat_groups=[(result.__dict__['structs'][0]['group'], result.__dict__['structs'][1]['group'])
                       for result in test_results]
                
               
                stat_dict['Groups'][g]=dict(zip(stat_groups, stat_results))    
                
            else:
                print('some')
                ax.set_ylim([-2, 10])

            ax.set_xlabel(g, fontsize=labelsize)
            

            #b.axes.set_title("Title",fontsize=50)
            
            #ax.set_ylim([-2, 6])
            
            ax.set_xticks([])
           
            ax.tick_params(labelsize=ticksize)
            ylabel=gen_y_label(dts[2])
            
            print(ylabel, 'ylabel')
            ax.set_ylabel(ylabel, fontsize=labelsize)
            
            if axisindex>0:
                 ax.set_ylabel(' ')
            
            
            #ax.get_legend().remove()
            #adjust_box_widths(fig, 0.9)
       
        
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.9, 1), fontsize=labelsize)
        
        plt.savefig(to+'/'+x+y+'bylabel'+str(annotit)+'.svg', bbox_inches='tight')
        
        print(stat_dict)
        
    with open(to+'\\'+plotype+y+"stat.pkl", "wb") as fp:
        pickle.dump(stat_dict, fp)


# In[ ]:





# In[48]:


def ecdf(a):
    x, counts = np.unique(a, return_counts=True)
    cusum = np.cumsum(counts)
    return x, cusum / cusum[-1]

def plot_ecdf(a):
    x, y = ecdf(a)
    x = np.insert(x, 0, x[0])
    y = np.insert(y, 0, 0.)
    plt.plot(x, y, drawstyle='steps-post')
    plt.grid(True)
    plt.savefig('ecdf.png')
    
    return (x, y)


# In[49]:


def plot_what(data, compare, param, select):
    """ compare: what to comare possible ['Group', 'Label', 'Age', 'Well iD, ....]
    param: param to compare
    select: subgroup in compare
    
    return data, x, y, hue
    """
    if len(select)>0:
        data=data[data[compare].isin(select)]
    x=compare   
    y=param
    
    
    return (data, x, y)


# In[ ]:





# In[50]:


def outlier_MCbox(df, param, MC=True):
    
    mc=medcouple(df[param].values)
    
    Q1=np.quantile(df[param].values, 0.25)
    Q3=np.quantile(df[param].values, 0.75)
    IQR=iqr(df[param].values)
    
    if MC==True:
        
        if mc>0:

            low=Q1-(1.5*IQR*e**(-3.5*mc))

            high=Q3+(1.5 * e** (4*mc) * IQR)

        elif mc<0:

            low=Q1-(1.5*IQR*e**(-4*mc))

            high=Q3+(1.5 * e** (3.5*mc) * IQR)

        else:

            low=Q1-(1.5*IQR)

            high=Q3+(1.5*IQR)
    else:
        
        low=Q1-(1.5*IQR)

        high=Q3+(1.5*IQR)
        
                      
    return low, high


# In[ ]:





# In[51]:


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


# In[52]:


def networkburstslider(wellspikes, start, end, slide, window, threshold_scaler, electrode_order):
    
    
    slidesxshape=(end-start)//slide
    
    windowshape=slidesxshape-(window//slide)


    yshape=12 ###n of the electrodes

    arr=np.zeros([slidesxshape+2, 12])


    timeaxis=np.arange(start, end, slide)[:-int(window//slide)+2]



    timeaxisends=timeaxis+window

    #print(timeaxis, timeaxisends, 'timeaxeses')



    for chlabel in range(yshape):

        spikes=wellspikes[wellspikes['Channel Label']==chlabel]['Timestamp [µs]'].values


        spikebinned=bin(spikes, start, end, slide)

        arr[:, int(chlabel)]=spikebinned

    dfroll=pd.DataFrame(arr)



    ###Rolling MEAN (other stat when mean)
    # for each channel individually
    arrNspikes=dfroll.rolling(window//slide, closed='both', axis=0).mean().values[int(window//slide)-1:, :]



    #print(arrNspikes.shape, window//slide)

    Nspikes=np.sum(arrNspikes, axis=1) #np.percentile(arrNspikes, 75, axis=1) # Selecet the 4th most active electrode

    #threshold=np.percentile(Nspikes, 75)

    threshold=np.std(Nspikes)*threshold_scaler

    #threshold=10

    #p=figure()

    #p.line(range(len(Nspikes)), Nspikes)

    #p.line(range(len(Nspikes)), np.ones(len(Nspikes))*threshold, color='black' )



    #plt.hlines(y=threshold, xmin=0, xmax=len(Nspikes), color='red')

    #show(p)


    #print(threshold, 'threshold', arr.shape, Nspikes.shape)



    networkburstindexes=np.where(Nspikes>=threshold)[0]+1




    # (np.diff(networkburstindexes)!=2



    eachnetworkburstindex=np.split(networkburstindexes, np.where((np.diff(networkburstindexes)!=1)

                                                                 & (np.diff(networkburstindexes)!=2) & (np.diff(networkburstindexes)!=3))
                                                                 [0]+1)


    channelsnetworkburstindex=[[np.where(arrNspikes>threshold)[1] for r in ar] for ar in eachnetworkburstindex]


    networkburstspikesNchannels=np.array([np.sum(arr[ar[0]:ar[-1]+1+window//slide, :], axis=0) for ar in eachnetworkburstindex])

    channelspertfnetworkburstindex=[[np.where(arrNspikes[r, :]>threshold/12)[0] for r in ar] for ar in eachnetworkburstindex]

    channelspernetworkburstindex=[np.unique(np.where(arrNspikes[ar, :]>threshold/12)[1])  for ar in eachnetworkburstindex]

    df=pd.DataFrame(networkburstspikesNchannels, columns=range(12))

    networksranges=[[timeaxis[ar[0]], timeaxisends[ar[-1]+1]] for ar in eachnetworkburstindex]




    networksspikes=[wellspikes[(wellspikes['Timestamp [µs]']>t[0]) 
                               & (wellspikes['Timestamp [µs]']<=t[1])].sort_values(by='Timestamp [µs]')['Timestamp [µs]'].values 
                    for t in networksranges]


    networkschannels=[wellspikes[(wellspikes['Timestamp [µs]']>t[0]) 
        & (wellspikes['Timestamp [µs]']<=t[1])].sort_values(by='Timestamp [µs]')['Channel Label'].values for t in networksranges]

    df['start']=[i[0] for i in networksranges]

    df['end']=[i[1] for i in networksranges]

    df['Duration']=[(i[1]-i[0])/1000 for i in networksranges]

    df['NetworkChannels']=channelspernetworkburstindex



    df['Number of Spikes']=[len(array) for array in networksspikes]


    df['MinISI']=[min(np.diff(sorted(array))) for array in networksspikes]

    df['MaxISI']=[max(np.diff(sorted(array))) for array in networksspikes]

    df['MeanISI']=[np.mean(np.diff(sorted(array))) for array in networksspikes]

    df['MedianISI']=[np.median(np.diff(sorted(array))) for array in networksspikes]

    df['MedianMeanISI']=df['MedianISI'].values/df['MeanISI'].values
    
    df['Network Spikes']=networksspikes
    
    df['Channels on Spikes']=networkschannels
    df['Network Ranges']=networksranges
    
    df['Detection Threshold']=threshold
    
    
    
    
    
    return df


# In[53]:


def networkburstslider2(wellspikes, start, end, slide, window, threshold_scaler, electrode_order):
    
    
    
    
    slidesxshape=(end-start)//slide
    
    windowshape=slidesxshape-(window//slide)


    yshape=12 ###n of the electrodes

    arr=np.zeros([slidesxshape+2, 12])


    timeaxis=np.arange(start, end, slide)[:-int(window//slide)+2]
    
    print(timeaxis.shape)



    timeaxisends=timeaxis+window

    #print(timeaxis, timeaxisends, 'timeaxeses')



    for chlabel in range(yshape):

        spikes=wellspikes[wellspikes['Channel Label']==chlabel]['Timestamp [µs]'].values
        #print(len(spikes), 'lenspikes')


        spikebinned=bin(spikes, start, end, slide)

        arr[:, int(chlabel)]=spikebinned
        
        
    #print(arr.shape, np.unique(arr), 'uniquarr')

    dfroll=pd.DataFrame(arr)



    ###Rolling MEAN (other stat when mean)
    # for each channel individually
    arrNspikes=dfroll.rolling(window//slide, closed='both', axis=0).mean().values[int(window//slide)-1:, :]
    
    #plt.plot(arrNspikes)
    
    #values[int(window//slide)-1:, :]



    

    Nspikes=np.sum(arrNspikes, axis=1) #np.percentile(arrNspikes, 75, axis=1) # Selecet the 4th most active electrode

    #threshold=np.percentile(Nspikes, 75)

    threshold=np.std(Nspikes)*threshold_scaler

    print(threshold, 'threshold')

    #p=figure()

    #p.line(range(len(Nspikes)), Nspikes)

    #p.line(range(len(Nspikes)), np.ones(len(Nspikes))*threshold, color='black' )



    #plt.hlines(y=threshold, xmin=0, xmax=len(Nspikes), color='red')

    #show(p)


    #print(threshold, 'threshold', arr.shape, Nspikes.shape)



    networkburstindexes=np.where(Nspikes>=threshold)[0]
    
    print(arrNspikes.shape, window//slide, Nspikes.shape)
    
    
        




    eachnetworkburstindex=np.split(networkburstindexes, np.where((np.diff(networkburstindexes)!=1)

                        & (np.diff(networkburstindexes)!=2) & (np.diff(networkburstindexes)!=3))[0]+1)
    
    
    #print(eachnetworkburstindex, 'eachnetworkburst')
    
   


    channelsnetworkburstindex=[[np.where(arrNspikes>threshold)[1] for r in ar] for ar in eachnetworkburstindex]


    networkburstspikesNchannels=np.array([np.sum(arr[ar[0]:ar[-1]+1+window//slide, :], axis=0)
                                          for ar in eachnetworkburstindex 
                                          if len(ar)>0])

    channelspertfnetworkburstindex=[[np.where(arrNspikes[r, :]>threshold/12)[0] for r in ar[:-1]]
                                    for ar in eachnetworkburstindex if len(ar)>0]

    channelspernetworkburstindex=[np.unique(np.where(arrNspikes[ar[:-1], :]>threshold/12)[1])  
                        for ar in eachnetworkburstindex if  len(ar)>0]

    df=pd.DataFrame(networkburstspikesNchannels, columns=range(12))

    networksranges=[[timeaxis[ar[0]], timeaxisends[ar[-1]-1]] for ar in eachnetworkburstindex if len(ar)>0]
    
    print('threshold', threshold)




    networksspikes=[wellspikes[(wellspikes['Timestamp [µs]']>t[0]) 
                               & (wellspikes['Timestamp [µs]']<=t[1])].sort_values(by='Timestamp [µs]')['Timestamp [µs]'].values 
                    for t in networksranges]
    
    #print([len(array) for array in networksspikes if len(array)==1], 'timeaxis, networkspikes' )


    networkschannels=[wellspikes[(wellspikes['Timestamp [µs]']>t[0]) 
        & (wellspikes['Timestamp [µs]']<=t[1])].sort_values(by='Timestamp [µs]')['Channel Label'].values for t in networksranges]

    df['start']=[i[0] for i in networksranges]

    df['end']=[i[1] for i in networksranges]

    df['Duration']=[(i[1]-i[0])/1000 for i in networksranges]

    df['NetworkChannels']=channelspernetworkburstindex



    df['Number of Spikes']=[len(array) for array in networksspikes]
    
    
    #print([len(array) for array in networksspikes], start, end, max(wellspikes['Timestamp [µs]'].values), min(wellspikes['Timestamp [µs]'].values))
    

    df['MinISI']=[min(np.diff(sorted(array))) if len(array)>1 else -1 for array in networksspikes ]

    df['MaxISI']=[max(np.diff(sorted(array))) if len(array)>1 else -1 for array in networksspikes]

    df['MeanISI']=[np.mean(np.diff(sorted(array))) if len(array)>1 else -1 for array in networksspikes]

    df['MedianISI']=[np.median(np.diff(sorted(array))) if len(array)>1 else -1 for array in networksspikes]

    df['MedianMeanISI']=df['MedianISI'].values/df['MeanISI'].values
    
    df['Network Spikes']=networksspikes
    
    df['Channels on Spikes']=networkschannels
    df['Network Ranges']=networksranges
    
    df['Detection Threshold']=threshold
    
    
    
    
    
    return networksspikes


# In[54]:


def rate_prop(dfall, percsval, param):
    
    """dfall is a grouped by object, percs is a list of percentiles of which proportions are caclulated"""
    
    
    rts=dfall[param].values
    
    
    df=pd.DataFrame()
    
   
    
    
    
    for pi in range(1, len(percsval)):
        
        rangeval=len(rts[(rts<percsval[pi]) & (rts>=percsval[pi-1])])
        
        
        
        
        df[str(int(percsval[pi-1]))+str(int(percsval[pi]))]=[rangeval/len(rts)]
        
    rangeval=len(rts[rts<percsval[0]])
    
    df[str(int(percsval[0]))+str(int(percsval[0]))]=[rangeval/len(rts)]
    
    rangeval=len(rts[rts>percsval[-1]])
    
    df[str(int(percsval[-1]))+str(int(percsval[-1]))]=[rangeval/len(rts)]
        
    
        
        
        
        
    return df


# In[55]:


def rate_prop(dfall, percsval, param):
    
    """dfall is a grouped by object, percs is a list of percentiles of which proportions are caclulated"""
    
    
    rts=dfall[param].values
    
    
    df=pd.DataFrame()
    
   
    
    
    
    for pi in range(1, len(percsval)):
        
        rangeval=len(rts[(rts<percsval[pi]) & (rts>=percsval[pi-1])])
        
        
        
        
        df[str(int(percsval[pi-1]))+str(int(percsval[pi]))]=[rangeval/len(rts)]
        
    rangeval=len(rts[rts<percsval[0]])
    
    df[str(int(percsval[0]))+str(int(percsval[0]))]=[rangeval/len(rts)]
    
    rangeval=len(rts[rts>percsval[-1]])
    
    df[str(int(percsval[-1]))+str(int(percsval[-1]))]=[rangeval/len(rts)]
        
    
        
        
        
        
    return df


# In[56]:


def load_meta(dpathexperiment):
    
    """dpathexperiment is a r string where your experimental spike files are located """
    
    stamp=np.load(dpathexperiment+'/stamp.npy', allow_pickle=True)[()]
    
    LED=np.load(dpathexperiment+'/LED_info.npy', allow_pickle=True)[()]
    
    Book=np.load(dpathexperiment+'/Lab_Book.npy', allow_pickle=True)[()]
    
    filenamegroup=[f for f in os.listdir(dpathexperiment) if 'layout' in f][0]
    
    Group_iN=np.load(dpathexperiment+'\\'+filenamegroup, allow_pickle=True)
    
    
    return stamp, LED, Book, Group_iN


# In[57]:


def weightrates(chrate):
    
    #rts[0]*(rts[1]/sum(rts[1])) was like this before, does it have any meaning
    
    rts=np.unique(chrate['spike rates'].values, return_counts=True)
    spikeuniqweighted=rts[0]*(rts[1]/sum(rts[1]))
    
    return spikeuniqweighted
    


# In[58]:


def return_bokeh_colormap(name):
    cm = plt.cm.get_cmap(name)
    colormap = [rgb_to_hex(tuple((np.array(cm(x))*255).astype('int'))) for x in range(0,cm.N)]
    return colormap
def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb[0:3]


# In[59]:


def relative_here(dburst, tosub, params):
    
    
    
    ss=dburst.groupby(['Experiment', 'Compound', 'Group'])['Label'].value_counts().min()
    ndburstdir=dburst.groupby(['Experiment', 'Compound', 'Group', 'Label']).apply(lambda x: 
                                                x.iloc[random.sample(range(0, len(x)), ss), :]).reset_index(drop=True)
    #params=['Duration', 'Spike Count', 'Max ISI', 'Min ISI', 'Mean ISI',
        #'Burstiness', 'FT Spikes',
        #'Surprise']
    
    
    ndburst=ndburstdir.groupby(['Experiment', 'Compound', 'Group']).apply(lambda x:
                                                                        relative_b(x, tosub, x['hueorder'].values[0],
                                                                    params, 'Group',  tosub, 'False').reset_index(drop=True))
        
    
 
        
        
    return (ndburst.reset_index(drop=True), ndburstdir)


# In[60]:


def recursive_items(dictionary):
    
    val=[]
   
    
    for key, value in sorted(dictionary.items()):
        
        
        if type(value) is dict:
            
            
            return recursive_items(value)
            
        else:
            
            for key, value in sorted(dictionary.items(), key=lambda item: item[1]):
            
                val.extend(value)
            
    return val


# In[61]:


def scale_df(df, cols):
    
    from sklearn.preprocessing import StandardScaler
    
    
    
    df[cols]= StandardScaler().fit_transform(df[cols])
    
    return df


# In[62]:


def load_hue_order(basepath):
    
    hueordp=[p[0] for p in os.walk(basepath)  for f in list(os.listdir(p[0])) if f=='hue_order.npy']
    
    hue_ord=np.load(hueordp[0]+'\\'+'hue_order.npy', allow_pickle=True)[()]
        
    return hue_ord


# In[63]:


def gen_y_label(y):
    
    if y.find('of')>-1:
        
        main=y[y.find('of')+2:]
        
        if main.find('Hz')>-1:
            
            res='Change'+' '+main[:main.find('(Hz)')]
        else:
            
            res='Change'+ ' '+main
    else:
        res=y
        
    return res


# In[64]:


def stat_testpub(data, x, hue, status):
    
    
    if len(np.unique(data[x].values).tolist())>1  and x!='Label' and len(np.unique(data[hue].values).tolist())>1:

        pair=Xhue(data, x, hue, intra=True, series=True)[1]
        #print(pair, "????")
        
    else:
        pair=Xhue(data, x, hue, intra=True, series=True)[0]

        #print(pair, "????new")
        
        
    if status=='dep':
        
        test1='Wilcoxon'
    else:
        test1='Mann-Whitney'
        
        
    if x=='Age':
        
        pair=Xhue(data, x, hue, intra=False, series=True)[0]
    
        
        
        
    return pair, test1


# In[65]:


def scatter_dosespub(dfs, experiment, Param, sortby, to, order, colordict, plotparam):
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
    
    
    print(order)
    
    combin=list(combinations(Labels, 2))
    
    combin=[comb[::-1]  if  np.where(order==comb[0])[0]>np.where(order==comb[1])[0] else comb for comb in combin]
    
    print(combin)
            
    
    
    combin= [comb for comb in combin if  np.where(order==comb[1])[0]-np.where(order==comb[0])[0]==1]
    
   
         
    combin=[combin[np.where(np.array(combin)[:, 0]==order[i])[0][0]] for i in range(len(order)-1)]
    
    print(combin, 'combin')
    
    
    
    fig, axs=plt.subplots(len(groups), len(combin), figsize=plotparam['figsize'],  sharey='col', sharex='col')
    
    if len(combin)==1:
        
        axs=axs.reshape(axs.shape[0], 1)
    
    for cindex, comb in enumerate(combin):
        
        #dfdoses=dfs[(dfs['Label']==comb[0]) | (dfs['Label']==comb[1])]
        
        #sns.scatterplot(data=dfdoses, x=comb[0], y=comb[1], hue='Compound')
        
        for index, comp in enumerate(comps):
            
            
            for indg, g in enumerate(groups):


                df=dfs[(dfs['Compound']==comp) & (dfs['Group']==g)]



                x=df[df['Label']==comb[0]].sort_values(by=sortby)

                colorlist=np.repeat(colordict[g], len(x))




                #markerswell=x['Well ID'].values.tolist()
                #markers=[markerlist[int(i)] for i in markerswell]

                #colormap=np.arange(0, len(x['Channel ID'].values), 1)



                x=x[Param].values
                y=df[df['Label']==comb[1]].sort_values(by=sortby)
                y=y[Param].values
                
                
                
                



                axs[indg, cindex].scatter(x, y, color=colordict[g], label=g)
                
                #sns.scatterplot(data=x, x=Param, y=, color=colordict[g], label=g)
                
                axs[indg, cindex].set_ylabel(comb[1],  fontsize=plotparam['labelsize'])
                #axs[indg, cindex].set_xlim([0, 20])
                #plt.xlim(0, max_value)
                #plt.ylim(0, max_value)





                model = LinearRegression().fit(x.reshape(-1, 1), y)#sm.OLS(y, sm.add_constant(x)).fit()
                y_hat=model.predict(x.reshape(-1, 1))

                r2=r2_score(y, y_hat)   

                #plt.plot(x, y_hat, color=colormap[index], label="r\u00b2 = {:.3f}".format(r2)+" "+comp)
                axs[indg, cindex].plot(x, x, color='black', label='Unity', linestyle='dotted')
                
                #axs[indg, cindex].set_xlim([-0.5, 50])


                #plt.annotate("r-squared = {:.3f}".format(r2), (0, 1), color=colormap[index])


                #plt.legend(fontsize=15)
                #plt.savefig(to+'/'+experiment+'Scatter'+comb[0]+comb[1]+comp+g+'.png')

                #influence = model.get_influence()
                #inf_sum = influence.summary_frame()


                #student_resid = influence.resid_studentized_external

       
        

        axs[indg, cindex].set_xlabel(comb[0], fontsize=plotparam['labelsize'])
        
        axs[0, 0].set_title(Param+experiment)
        handles=[mpatches.Patch(color=v, label=k) for k, v in colordict.items()]
            
        fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.9, 1), fontsize=plotparam['labelsize'])
       
        plt.savefig(to+'/'+experiment+Param+'Scatterxlim'+'.png', bbox_inches='tight')
        
    return Labels, combin


# In[66]:


def meltgroup(dburstmedians, burstpcols, group1, group2):
    dburstmediansW=dburstmedians.groupby(['Base', 'Group', 'Well ID', 'Week']).apply(lambda x: x[burstpcols].mean()).reset_index()
    dburstmediansWgroups=dburstmediansW.groupby('Base').filter(lambda x: ((group1 in x['Group'].values)==True)
                                                                     & ((group2 in x['Group'].values)==True))


    dburstmediansWgroupscaled=dburstmediansWgroups.groupby(['Week', 'Base']).apply(lambda x: 
                                                                              scale_df(x[(x['Group']==group1) | 
                        (x['Group']==group2)], burstpcols)).reset_index(drop=True)  


    dburstmediansWgroupscaledmelted=pd.melt(dburstmediansWgroupscaled, id_vars =['Base', 'Group', 'Well ID', 'Week'],
                                            value_vars =burstpcols, value_name='Z scores') 
    
    return  dburstmediansWgroupscaled, dburstmediansWgroupscaledmelted




def plot_lineplot(to, datasetmelted, I1, I2, palette, hue_order,  estimator='mean', 
                  errorbar=('se', 2), rangeval=2.5, figsize=(6, 4), style='Week'):



    plt.figure(figsize=figsize)

    sns.set(font_scale=1.1)
    sns.lineplot(data=datasetmelted[(datasetmelted['Z scores']>-2.5) & (datasetmelted['Z scores']<2.5)], 
                 x='variable', y='Z scores', hue='Group', palette=palette, hue_order=hue_order, style=style, errorbar=errorbar, estimator=estimator)
    
    plt.title(I1+I2)


    plt.savefig(to+'/'+I1+I2+'Line11823'+ errorbar[0]+estimator + '.png')
    
    plt.show()
    
def boxbar(df, AgeParam, param, compare_order, colour_dict):

    signpvalues=df.groupby(AgeParam).apply(lambda x: statistical_test(x, 'Group', param, 'Ch', stats.mannwhitneyu, False,  [])[1])
    signpairs=df.groupby(AgeParam).apply(lambda x: statistical_test(x, 'Group', param, 'Ch', stats.mannwhitneyu, False,  [])[2])
    dit=list(signpairs.to_dict().items())

    pvaldict=list(signpvalues.to_dict().items())
    pairs=[((pair[0],  p[0]), (pair[0], p[1]))  for pair in dit for p in pair[1]]
    pvalues=[pvaldict[i][1][p] for i in range(len(dit)) for p in dit[i][1]]

    plt.figure(figsize=(6, 4))

    ax=sns.barplot(data=df, x='Age', y=param, hue='Group', hue_order=compare_order, palette=colour_dict.values())


    annot = Annotator(ax, pairs, data=df, x='Age',  y=param, hue='Group',
                   hue_order=compare_order, palette=colour_dict.values(),  boxprops=dict(alpha=.7))#, order=order, hue_order=hue_order)


    #annot.configure(test=test1, text_format='star', loc='outside', comparisons_correction="bonferroni")
    annot.configure(loc='outside')


    annot.set_pvalues(pvalues)


    plt.legend(bbox_to_anchor=(1, 1))


    annot.annotate()


    #plt.title(df['Base'].unique()[0])

    plt.show()


# # Spike Labeling Functions

# In[4]:


def Spike_Source_Details(mergeddf, column):
    
    '''column is bursting or network bursting, where 3 means found in both spikes and col and 2 only spikes
    merfeddf has Bursting or Network Bursting columns'''
    
    mergeddf=mergeddf.sort_values(by='Timestamp [µs]')
    
    Bursting=mergeddf[column].values.astype('int')
    
    spiketimes=mergeddf['Timestamp [µs]'].values

    Diffbursting=Bursting[1:]-Bursting[:-1]
    
    ISI=np.log(spiketimes[1:]-spiketimes[:-1])
    
    icb=np.where((Diffbursting==0) & (Bursting[1:]==3))[0]
    ics=np.where((Diffbursting==0) & (Bursting[1:]==2))[0]

    its=np.where((Diffbursting==-1))[0]
    itb=np.where((Diffbursting==1))[0]

    arr=np.zeros(len(Diffbursting)).astype('str')
    arr[icb]='In Burst'
    arr[ics]='Irregular'
    arr[its]='To irregular'
    arr[itb]='To Burst'
    
    arr=[arr[0]]+arr.tolist()
    
    mergeddf[column+'Details']=arr
    
    mergeddf['ISI']=[ISI[0]]+ISI.tolist()
    
    return mergeddf


# # Group of function to check LED stimulation condition and params if any. 

# In[67]:


def led_param(file, led, df, stamp, listamp, label_index, dlabel):

    '''led is LED dict, df is dataset grouped by exp and dose , stamp is times of given exp at dose
    listamp is user input for possible amps in experiment, as no amplitude info was found  amp should be in file name'''
    
    pos=file.find('h5')+2
    
    
    
    key1=file[:pos]
    key2=dlabel
    
    dose_start, dose_end=segment(df, stamp, key2)
    
    
    
    start, end, leds, led_wells, light_label, light_dur=dose_end, dose_end, [], [], '', 0
    
    ##label_index defines index for specifc light 
    times=np.array(led[key1]['Time'][label_index])
    
    ##choose led times that belong to a given dose label 
    
    if len(times[(times>dose_start) & (times<dose_end)])>0:
        
        ##leds within the dose 
        
        leds=times[(times>=dose_start) & (times<=dose_end)]
        
        start=leds[0]
        
        
        end=leds[-1]
    
    
    
    #print(start, end)
    
    
        LED_labels=led[key1]['Wells'][label_index]

       

        led_wells=label_to_id(LED_labels)[0]
        light_label=led[key1]['Label'][label_index]

        light_label=[string[-6:].lower() for string in light_label]

        light_dur=led[key1]['Duration'][0][0] #in microseconds as everything else 



        light_freq=np.mean(1000000/np.diff(leds))



        if len(listamp)>0:



                amps=find_amp_weird(df['Experiment'].unique()[0],  listamp)


        else:
                amps=14

    
    return (start, end, leds, led_wells, light_label, light_dur) 





# In[68]:


def label_to_id(led_label):
    
    #the led and stimulation are given by well labels, we work on ids. 
    
    WellID=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellID.csv.npy")
    Welllabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\WellLabel.csv.npy")
    Channellabel=np.load(r"C:\Users\MEA_PC\Desktop\AH\Python\Metadata\ChannelLabel.csv.npy")
    
    
    led_id=[]
    
    led_label=led_label.split(',')
    l_condition=np.repeat('off', 24)
    for label in led_label:
        
        
        
        i=int(WellID[np.where(Welllabel==label)[0][0]])
        l_condition[i]='On'
        led_id.append(i)
        
    return led_id, l_condition


# In[69]:


def any_led(file, led, stamp, df, dlabel):
    
    anyarg=False
    
    led_arg=False
    
    pos=file.find('h5')+2
    
    
    key1=file[:pos]
    
    key2=dlabel
    
    dose_start, dose_end=segment(df, stamp, key2)
    
    if key1 in list(led.keys()):
        
        if len(list(led[key1].keys()))>0:
            
            if len(led[key1]['Time'])>0:
                
                anyarg=True
                
    return anyarg


# In[2]:


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
    
    


# In[71]:


def led_seg(file, led, df, stamp, dlabel, label_index):
    
   
    
    pos=file.find('h5')+2
    
    
    
    key1=file[:pos]
    key2=dlabel
    
    dose_start, dose_end=segment(df, stamp, key2)
    
    times=np.array(led[key1]['Time'][label_index])
    
    if len(times[(times>dose_start) & (times<dose_end)])>0:
        
        start=times[(times>dose_start) & (times<dose_end)][0]
        end=times[(times>dose_start) & (times<dose_end)][-1]
    
    return (start, end)


# In[72]:


def check_led(file, led, stamp, df, label_index, dlabel):
    
    
    led_arg=False
    
    
    
    pos=file.find('h5')+2
    
    
    key1=file[:pos]
    
    key2=dlabel
    
    dose_start, dose_end=segment(df, stamp, key2)
    
    if key1 in list(led.keys()):
        
        if len(list(led[key1].keys()))>0:
            
            if len(led[key1]['Time'])>=label_index+1:
                
                times=np.array(led[key1]['Time'][label_index])
                
                if len(times[(times>dose_start) & (times<dose_end)])>0:
                    
                    led_arg=True
                    
    return led_arg
                   


# # Chanel selection, distributions and  so on. 

# In[73]:


def clean_bad_channels(df, param):
    
    '''param is BaselineM'''

    
    high=outlier_MCbox(df, param)
    
    channels_good=np.unique(df[df[param]<high]['Channel ID'].values)
    
    return channels_good


# In[ ]:





# 
# # LED and PSTH functions 
# 

# In[ ]:





# In[5]:


def relativeISI_merge(data, col, relativekey1, relativekey2, param):
    
    
    
    ###if mergeon not really needed to merge

    mergeoncols=['Channel ID', 'Well ID', 'Well Label', 'Channel Label', 
           'Dose Label', 'Compound ID', 'Group', 'Time_Miliseconds', 'Induction', 'LED Start']

    pre=data[data[col]==relativekey1].sort_values(['Channel ID', 'LED Start'])[mergeoncols+[param]]
    post=data[data[col]==relativekey2].sort_values(['Channel ID', 'LED Start'])[mergeoncols+[param]]


    dictLEDstartmap=dict(zip(sorted(pre['LED Start'].unique()), sorted(post['LED Start'].unique())))

    pre['LED Start']=pre['LED Start'].map(dictLEDstartmap)

    merged=pd.merge(pre, post, on=mergeoncols, how='outer', suffixes=['pre', 'post'])



    merged.loc[:, param+'_pct']=(merged['1/ISIpost']-merged['1/ISIpre'])/(merged['1/ISIpre']+0.00001)
    
    
    return merged


# # Calculation of rates, including instanteous, Average, binned Rates and Statistics
# 

# In[1]:


def calc_micro_rate_for_LED(df, LED, stamp, winsize=500, window=5000, bfparam=300, afparam=500):
    
    
    '''calculates led evoked spiking insta rates on grouped by dataset, 
    for 0.5ms bins and for 5ms, df is a grouped by dataset'''
    
    
    
     
    df_collect2ms=pd.DataFrame()
    
    df_collect5ms=pd.DataFrame()
    
    if len(df['Dose Label'].unique())==1:
    
    
    
        dlabel=df['Dose Label'].unique()[0]
        
        file=df['File'].unique()[0]


        wellid=df['Well ID'].unique()[0]


        

        start, end=segment(df, stamp, dlabel, Duration=None)

        ##micro rates from 0 to the end
        
        
    


        lightlabels=led_label(file, LED, df, stamp, dlabel)
        
        
        if len(lightlabels)>0:



            for index, name in enumerate(lightlabels):


                ledarg=check_led(file, LED, stamp, df, index, dlabel)

                if ledarg==True:

                    print(ledarg)

                    #500 us binned data

                    arr=micro_rate(df, end, winsize)


                    ledstart, ledend, leds, led_wells, light_label, light_dur=led_param(file, LED, df, 
                                                                                        stamp, [], index, dlabel)
                   
                    if len(leds)>1:

                        interled=np.diff(leds)[0]
                    else:
                        interled=1000000
                   

                    print(wellid, led_wells, len(leds))

                    if wellid in led_wells:

                        print('well id ')


                        for ld in leds:

                            print(ld,  'led', start, end, 'start, end')

                            ###return is correct is returns the correct lde 

                            df5ms, df2ms=led_micro_rate(arr, ld, interled, winsize, window, bfparam, afparam)


                            df2ms['LighWaveform']=name
                            df2ms['LED Start']=ld

                            df2ms['Light Frequency']=1000000/interled


                            df5ms['LighWaveform']=name
                            df5ms['LED Start']=ld

                            df5ms['Light Frequency']=1000000/interled




                            df_collect2ms=pd.concat([df_collect2ms, df2ms])
                            df_collect5ms=pd.concat([df_collect5ms, df5ms])

            else:
                print('no led, empty dataset')

                
    return (df_collect2ms, df_collect5ms)


# In[4]:


def Calculate_Spontant_Micro_rates(df, stamp, led, winsize=500, window=5000):
    
    dlabel=df['Dose Label'].unique()[0]

    file=df['File'].unique()[0]


    wellid=df['Well ID'].unique()[0]

    start, end=segment(df, stamp, dlabel, Duration=None)
    
    ledarg=check_led(file, led, stamp, df, 0, dlabel)
    
    arr=micro_rate(df, end, winsize)
    
    if ledarg==False:
        
        
        
        smoothed_df=smooth_micro_rate(arr, start, end, winsize=500, window=5000)
        
        
    else:
        
        ledstart, ledend, leds, led_wells, light_label, light_dur=led_param(file, led, df, 
                                                                                        stamp, [], 0, dlabel)
        
        
        
        smoothed_df=smooth_micro_rate(arr, start, ledstart, winsize=500, window=5000)
        
        
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
               
               
         


# In[2]:


def led_micro_rate(ar, led, interled, winsize, window, bfparam, afparam):
    
    
    '''In order to avoid inprecision in LED time led arrays are smoothed separately if windoe is 5000 snd 500
    last 10 vslues are left, baseline 100ms is taken 
    
    takes led and laready 500us microrate array and smoothes 
    
    around a windows arround led stimulus, considering baseline and post induction window for 1000 ms '''
 
    ##start is including the baseline and can add plotting time.
    
   

    ms2df=pd.DataFrame({'1/ISI': [], 'Time': [], 'Time_from_start':[],  
                            'Time_to_plot':[]})
    
    ##in addition to that i can have 100ms beofre the interval 
    
    beforeled=bfparam*1000 #300 ms
    afterled=(afparam*1000)# 500 ms 
   
    start=led-beforeled
    
    end=led+afterled
    
    
    
    if start < 0:
        
        print('no baseline 100ms')
        
    else:
        
        print('with baseline 100ms')
        
        
        timeax=np.linspace(start, end, len(ar[int(start//winsize):int(end//winsize)]))-start
        
       

        timeaxtrue=timeax+start
        
        ###plotting var that start witn negative value. 
        
        timeaxplotting=(timeax-beforeled)//1000
        
        
        #500us spikes in range of beforeled-postled
        
        arinrange=ar[int(start//winsize):int(end//winsize)]
        
        


       
        
        dfroll=pd.DataFrame(arinrange,
                            index=pd.TimedeltaIndex(timeax, unit='us'), columns=['1/ISI'])
        
        ###rolling per 2ms but subsampling for 5ms , but good param needs to be checked

        arr=dfroll.rolling('2000us', closed='both', axis=0).mean().resample(str(window)+'us').mean()
        
        
        ###considering that starts from baseline, plot as negative

        arr['Time_to_plot']=arr.index.astype('timedelta64[ms]').values-(beforeled//1000)
        
       # print(arr.shape,   dfroll.rolling('2000us', closed='both', axis=0).mean().shape, str(window)+'us', 'shapes before, after resample')
        
        
        #time_reference=np.arange(beforeled, afterled, window)[:arr.shape[0]]
        
        #arr['Time_Reference']=time_reference
        
        

        
        arr['Time_from_start']=arr.index.astype('timedelta64[us]').values+start
        
        arr['Time']=arr.index.astype('timedelta64[us]').values
        

        ##Rolling MEAN (other stat when mean)
        #for each channel individually


       
        ###m2df, 500us, 

        ms2df=pd.DataFrame({'1/ISI': arinrange, 'Time': timeax, 'Time_from_start':timeaxtrue, 
                            
                            'Time_to_plot':timeaxplotting})
        
       
        
      

    
    return arr, ms2df 


# In[82]:


def calcspontanmfr(df, led, stamp, Param, post=False):
    
    
    ##it is different when func runs in a groupby object or on df. Df has 1 unique object
    
    mfr=-1
    
    
    if len(df['Dose Label'].unique())==1:
        
        file=df['File'].unique()[0]
        
        
    
        Label=df['Dose Label'].unique()[0]


        #when to start after induction 
        postinduction=1000000

        start, end=segment(df, stamp, Label, Duration=None)

        ##checking leds

        anyled=any_led(file, led, stamp, df, Label)

        if anyled==True:

            label_index=led_label(file, led, df, stamp, Label)

            if len(label_index)>=1:

                led_arg=check_led(file, led, stamp, df, 0, Label)

                if led_arg==True:

                    leds, lede=led_seg(file, led, df, stamp, Label, 0)
                    ##if post stimulation spontant, starts after
                    if post==True:

                        start=lede+postinduction
                    #else  whole calculation ends before LED starts     
                    else:

                        end=leds


        duration=(end-start)/1000000 # duration of recording in sec 

        spikes=df[Param].values

        spikes=spikes[(spikes>start) & (spikes<end)]

        mfr=len(spikes)/duration


    
    return mfr
        


# In[ ]:





# # Plotting objects

# In[83]:


def Set_Seaborn_Styles():
    
   
    # This sets reasonable defaults for font size for
    # a figure that will go in a paper
    sns.set_context("paper")
    
    # Set the font to be serif, rather than sans
    sns.set(font='serif')
    
    # Make the background white, and specify the
    # specific font family
    sns.set_style("white", {
        "font.family": "serif",
        "font.serif": ["Times", "Palatino", "serif"]
    })
    
  


# In[ ]:





# In[4]:


def plot_facegrid(data, snsplottype, plotdict):
    # Specify that I want each subplot to correspond to
    # a different robot type
    Set_Seaborn_Styles()
    
    
    with sns.plotting_context(rc={"legend.fontsize":12}):
    
        g = sns.relplot(data=data, **plotdict)

        # Create the bar plot on each subplot

        axes = np.array(g.axes.flat)
        
        for a in axes:
            
            a.axvline(x=0)

            a.set_xlabel('Time(ms)')
    
       
        
    return plt.gcf(), axes


# In[3]:


def plot_facegrid1(data, snsplottype, x, y, colgrid, col_order,  sharey, sharex, hue,  hue_order):
    # Specify that I want each subplot to correspond to
    # a different robot type
    Set_Seaborn_Styles()
    
    
    with sns.plotting_context(rc={"legend.fontsize":12}):
    
        g = sns.relplot(data=data,
            x=x, y=y, hue=hue, hue_order=hue_order, kind='line',  col=colgrid,
            col_order=col_order, errorbar='ci')

        # Create the bar plot on each subplot

        #axes = np.array(g.axes.flat)
       
    
       
        
    return plt.gcf(), axes


# # Plotting multi-feature spider 
# 

# In[1]:


def normalize(value, param_range):
    min_val, max_val = param_range
    return (value - min_val) / (max_val - min_val) * 2 - 1 
def round_angle(angle, base=60):
    return int(base * round(angle / base))


def radar_plot(data, Params, my_palette, group_order,  group_col='Group'):
    
    
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

    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
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
    ax.set_ylim(-1.2, 1.2)
    #ax.set_ylim(-1, 1)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

   


# #  mapping raw to spikes

# In[ ]:





# In[6]:


def tk_get_raw(self, channelid, argumentsnew={}, durationofpart=5, methodarg=(np.median, [])):
    
        
        keys=['fs','cutoff', 'order','dead_time','threshold', 'negative', 'positive', 'mindur']
        
        arguments=self.spikeparams_set['spike detection kwargs']
        
        argumentsmethod={k: arguments[k] for k in keys if k not in ['order', 'cutoff']}
        
        argumentsmethod['durationofpart']=durationofpart
        
        for key, value in argumentsnew.items():
            
            argumentsmethod[key]=value
            
            
        folder_path=self.hdf5path
        
        file=self.rawfile
        
        dlabel=self.wellspikes['Dose Label'].unique()[0]
        
        print('Dose label raw', dlabel)
    
        Recordingfile=os.path.join(folder_path, file)
        
        print(Recordingfile, 'Recordign File')
    
        recording=McsPy.McsData.RawData(Recordingfile)
        
        
    
        index=int(recording.recordings[0].analog_streams[0].channel_infos[channelid].row_index)
        
        #print(index, recording.recordings[0].analog_streams[0].channel_infos[channelid].row_index )
        
        full_times=self.stamp[file][str([dlabel])]
        
        dose_st=full_times['start']
        dose_end=full_times['stop']
        
        starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
        stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]
        esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
        
        dose_start=starts[np.where(esim==dose_st)[0]][0]
        
        dose_stop=stops[np.where(esim==dose_st)[0]][0]
    
                  
        v=recording.recordings[0].analog_streams[0].channel_data[index, dose_start:dose_stop]
        
        fv=highpass(v, arguments['cutoff'], arguments['order'], arguments['fs'])
        
        spikes, theshold=detect_single(fv, methodarg, **argumentsmethod, start=dose_start, detection=False)
        
        ##dose_st must be +to spikes*50, as it is already dose_start*50
        
        return  (spikes, theshold, fv, dose_st, dose_end)
    
    
def tk_get_raw_traces(self, channelid, fv, bursts, dose_st,  threshold, methodstring):



    rawsxs=[]
    rawsys=[]


    ##500ms plot

    p1=figure(width=1000, height=200, title=str(channelid)+' '+str(threshold)+' '+methodstring)

    for key, values in bursts.items():

        for index, row in values.iterrows():

            spk=row['Spikes']


            times=[int(time_reverse(x, dose_st)) for x in spk]

            start=times[0]-100000

            if start<=0:

                start=0


            end=times[-1]+100000



            xval=np.arange(start, end, 1)

            rawsxs.append(xval.tolist())

            yval=fv[int(start):int(end)]

            rawsys.extend(yval.tolist())

    intrawsxs=[np.array(sub)-sub[0] for sub in rawsxs]  
    newrawsxspre=[np.array(sub)+intrawsxs[i-1][-1]+10000 if i!=0 else np.array(sub)
               for i, sub in enumerate(intrawsxs)]


    newrawsxs=[i for item in newrawsxspre for i in item]


    p1.line(newrawsxs, rawsys)

    hline1 = Span(location=threshold, dimension='width', line_color='red', line_width=2, line_dash='dotted')

    hline2 = Span(location=threshold*(-1), dimension='width', line_color='red', line_width=2, line_dash='dotted')

    p1.add_layout(hline1)

    p1.add_layout(hline2)


    self.rawplot_list.append(p1)


# # Selecting the data , filtering 

# In[2]:


###choose columns by abolsute median values and standard deviation, df must be changes dataset

def sort_describe_by_abs_median_and_std(df, top_n=10):
    # Get the describe output
    describe_df = df.describe()
    
    # Sorting by max absolute mean value and min std
    sorted_columns = describe_df.loc[['50%', 'std']].T
    sorted_columns['abs_median'] = sorted_columns['50%'].abs()

    # Sort first by absolute mean (descending) and then by std (ascending)
    sorted_columns = sorted_columns.sort_values(by=['abs_median', 'std'], ascending=[False, True])

    # Return the top N column names
    return sorted_columns.head(top_n).index.tolist()


# In[26]:


def condition_choose_query(df, condition_dict):
    
    query_str = ' & '.join([f'`{key}` == {repr(value)}' for key, value in condition_dict.items()])
    
    dfplot= df.query(query_str)
    
    return dfplot


# In[ ]:





# In[ ]:





# In[ ]:





# In[45]:


def clean_low_channels(df, param, threshold, selection):
    
    '''param is mfr, nbr, br, if network burst slection is well id'''

    
    
    
    selection_good=np.unique(df[df[param]>threshold][selection].values)
    
    return selection_good


# In[17]:


def channel_selecting(df, param, paramstat):
    
    """Using quantile values for all channels 
    for all well per day to select channels (""sensory that detects""), in general param should be 'mean rates'"""
    
   
    for q, name in zip([paramstat], ['S'+str(paramstat)]):
        
        channel_select=np.zeros([len(df)])
        
        indexes=np.where(df[param].values>q)[0]
        
        channel_select[indexes]=1
        
        df[name]=channel_select
        
    return df


# In[32]:


def ISIStat(df, timestampcol, param):
    
   
    spikes=df[timestampcol].values
    
    isi=np.diff(spikes)
    
    if len(isi)>1:

        df[param]=isi.tolist()+[isi[-1]]



        new_df=param_DistStat(df, param)

        return new_df


# In[54]:


def active_channels(df, sel_param, by):
    
    """function to calulate active subjct( channel or well) per group, per age, per well defiend by by"""
    
    
    
    df['Active Channels'+param]=len(df[df[sel_param]==1])
    
    
    return df


# In[61]:


def weightrate(chrate, param):
    
    '''chrates is the dataset with rate 
    param is rate param '''
    
    #rts[0]*(rts[1]/sum(rts[1])) was like this before, does it have any meaning
    
    rts=np.unique(chrate[param].values, return_counts=True)
    spikeuniqweighted=rts[0]*(rts[1]/sum(rts[1]))
    
    return spikeuniqweighted


# In[62]:


def WeightMFR(df):
    
    """WeigHt, df is already grouped by well"""
    
    
    df=df.sort_values(by=['Channel ID'], ascending=True)
    
    
    channels=sorted(df['Channel ID'].unique())
    
    meanrates=df['Mean of rates'].values
    
    summfr=sum(meanrates)
    
    wmfrs=[(mfr/summfr) for mfr in meanrates]
    
    df['WeightchMFR']=wmfrs
    
    return df
    


# In[66]:


def WMFR(df, param):
    
    """WEIGHTED MEAN FIRING CALCULATON, df is already grouped by well"""
    
    
    df=df.sort_values(by=['Channel ID'], ascending=True)
    
    
    channels=sorted(df['Channel ID'].unique())
    
    meanrates=df[param].values
    
    summfr=sum(meanrates)
    
    wmfrs=summfr/len(channels)
    
    
    
    return wmfrs
    


# In[2]:


def modify_b(b):
    
    if isinstance(b, str)==True:
    
        b=b[1:-1].split(' ')

        b=np.array([ch for ch in b if ch[:1].isnumeric()==True]).astype('float')
    
    
    return b


# In[2]:


def find_peak_slopev(b):
    
    
    
    b=modify_b(b)
    
    
    st=b[0]
    
    isis=np.diff(b)
    
    minisiind=np.argmin(isis)+1
    
    sloperise=(1/(isis/1000000))/((b[minisiind]-b[0])/1000000)
    
    slopefall=(1/(isis/1000000))/((b[-1]-b[minisiind])/1000000)
    
    
    
    return (sloperise, slopefall)


# # Network codes

# In[59]:


def networkold(wellb):
    
    channels=wellb['True Label'].values
    bursts=wellb['Timestamp [µs]'].values

    netwrokschs=[]

    netwroksbursts=[]

    indexall=[]
    
    

    for burst in bursts:



        arrdiff=np.array(bursts)-burst

        indexes=np.where(np.abs(arrdiff/1000) <100)[0]

        indexes=[ind for ind in indexes if ind not in indexall]

        indexall.extend(indexes)

        network=channels[indexes]

        netwroksburst=bursts[indexes]



        netwrokschs.append(network)

        netwroksbursts.append(netwroksburst)
        
    wellb['Network Bursts']=netwroksbursts
    
    wellb['Network Channels']=netwrokschs
    
    return wellb
 

    
    


# In[ ]:





# In[1]:


def MEA_toolboxnetworks():
    
# Check SCBflag for synchronization burst behavior
    if SCBflag == 0:
        # Efficiently collect all start times from bursts in burst7
        starttimesbursts2 = np.hstack([burst7[i][j][0] for i in range(len(burst7)) for j in range(len(burst7[i]))])

        # Unique start times of detected bursts
        uniquestarttimesbursts = np.unique(starttimesbursts2)

        # Proceed with further analysis on uniquestarttimesbursts...

    elif SCBflag == 1:
        # Clean up potentialnnbursts using mean and standard deviation thresholds
        potentialnnburstsmean = [
            np.nanmedian([np.nanmedian(burst) for burst in potentialnnbursts[i]]) + 
            2.5 * np.nanstd([np.nanstd(burst) for burst in potentialnnbursts[i]])
            for i in range(len(potentialnnbursts))
        ]

        # Eliminate bursts that are too far from the median
        for i in range(len(potentialnnbursts)):
            potentialnnbursts[i] = [
                burst for burst in potentialnnbursts[i] 
                if np.nanmedian(burst) <= potentialnnburstsmean[i]
            ]

        # Creating vectors for faster analysis (flatten burst lists)
        potentialnnburstsvector = [np.hstack(burst) for burst in potentialnnbursts]
        potentialnnburstsvector_min_max = np.array([[np.min(burst), np.max(burst)] for burst in potentialnnburstsvector])

        # Initialize lists for merging bursts
        tobemergedbursts = [[] for _ in range(len(potentialnnbursts))]
        mergeburstchannel = [[] for _ in range(len(potentialnnbursts))]
        mergeburstnumber = [[] for _ in range(len(potentialnnbursts))]

        # Vectorized search for bursts within synchronized windows
        for i, (start, end) in enumerate(potentialnnburstsvector_min_max):
            for k, channel_bursts in enumerate(burst6):
                for l, burst in enumerate(channel_bursts):
                    if np.any((burst > start) & (burst < end)):  # Vectorized comparison
                        tobemergedbursts[i].append(burst)
                        mergeburstchannel[i].append(k)
                        mergeburstnumber[i].append(l)

        # Clean up empty cells in a single pass
        tobemergedbursts = [b for b in tobemergedbursts if b]
        mergeburstchannel = [c for c in mergeburstchannel if c]
        mergeburstnumber = [n for n in mergeburstnumber if n]

        # Remove bursts with less than minimum channel participation
        min_participation = len([ch for ch in M2 if ch]) * fs['networkburstdetection']['minimumchannelparticipation']
        valid_indices = [i for i, channels in enumerate(mergeburstchannel) if len(channels) >= min_participation]

        # Apply the valid indices to filter bursts
        tobemergedbursts = [tobemergedbursts[i] for i in valid_indices]
        mergeburstchannel = [mergeburstchannel[i] for i in valid_indices]
        mergeburstnumber = [mergeburstnumber[i] for i in valid_indices]

        # Find and merge overlapping bursts
        findoverlappingnnbursts = np.array([[np.min(b), np.max(b)] for b in tobemergedbursts])

        # Merge overlapping bursts efficiently
        for i, burst_range in enumerate(findoverlappingnnbursts):
            if burst_range.size == 0:
                continue
            for k in range(i + 1, len(findoverlappingnnbursts)):
                if findoverlappingnnbursts[k].size == 0:
                    continue
                if (findoverlappingnnbursts[k][0] < burst_range[1] + 1) and (findoverlappingnnbursts[k][1] > burst_range[0] - 1):
                    # Merge bursts if overlap exists
                    findoverlappingnnbursts[i][1] = max(burst_range[1], findoverlappingnnbursts[k][1])
                    findoverlappingnnbursts[k] = np.array([])

        # Final filtering based on overlapping bursts
        valid_indices = [i for i, overlap in enumerate(findoverlappingnnbursts) if overlap.size > 0]
        tobemergedbursts = [tobemergedbursts[i] for i in valid_indices]
        mergeburstchannel = [mergeburstchannel[i] for i in valid_indices]
        mergeburstnumber = [mergeburstnumber[i] for i in valid_indices]
        findoverlappingnnbursts = [findoverlappingnnbursts[i] for i in valid_indices]

      
            
        
        return findoverlappingnnbursts


# In[51]:


def network(wellb, dtscale):
    
    '''Detects with 100ms threshold later can be filtered by sttc values  by multivariantly
    dtscale is 2ms more commonly'''
    
    
    wellb=wellb.reset_index()
    
    channels=wellb['True Label'].values
    bursts=wellb['Timestamp [µs]'].values
    
    spikes=wellb['Spikes'].values.tolist()

    netwrokschs=[]

    netwroksbursts=[]
    
    netwroksspikes=[]

    indexall=[]
    
    

    for burst in bursts:



        arrdiff=np.array(bursts)-burst

        indexes=np.where(np.abs(arrdiff/1000) <100)[0]

        indexesnew=[ind for ind in indexes if ind not in indexall]
        
        indexold=[ind for ind in indexes if ind in indexall]
        
        
        network=channels[indexesnew]

        netwroksburst=bursts[indexesnew]
        
        netwroksspike=[spikes[i] for i in indexesnew]#spikes[indexesnew]
        
        
        
        if len(indexold)>0:
            
            oldnetwork=bursts[indexold]
            tonetwroksbursts=[i for i in netwroksbursts if np.array(oldnetwork).any() in i]
            if len(tonetwroksbursts)>0:
                
                netwroksbursts=[i.extend(netwroksburst) for i in netwroksbursts if np.array(oldnetwork).any() in i]
                netwrokschs=[netwrokschs[i].extend(network) for i in  
                                range(len(netwroksbursts)) if np.array(oldnetwork).any() in netwroksbursts[i]]
                
                netwroksspikes=[netwrokspikes[i].extend(netwroksspike) for i in  
                                range(len(netwroksbursts)) if np.array(oldnetwork).any() in netwroksbursts[i]]
                
            else:
                netwroksbursts.append(netwroksburst)
                netwrokschs.append(network)
                netwroksspikes.append(netwroksspike)
        else:
            netwroksbursts.append(netwroksburst)
            netwrokschs.append(network)
            netwroksspikes.append(netwroksspike)
            
            
        indexall.extend(indexesnew)
        
        
    selectbursts=[ind for ind in range(len(netwrokschs)) if  len(np.unique(netwrokschs[ind]))>3]
    
    netwroksbursts=[netwroksbursts[ind] for ind in range(len(netwroksbursts)) if ind in selectbursts]
    
    netwrokschs=[netwrokschs[ind] for ind in range(len(netwrokschs)) if ind in selectbursts]
    
    netwroksspikes=[netwroksspikes[ind] for ind in range(len(netwroksspikes)) if ind in selectbursts]
    
        
    #netwroksbursts=[b for b in netwroksbursts if len(b)>3]
    
    #netwrokschs=[b for b in netwrokschs if len(b)>3]
    #netwroksspikes=[b for b in netwroksspikes if len(b)>3]
    
    allnetworkspikes=[np.unique(sorted([item for i in b for item in i ])).tolist() for b in  netwroksspikes]
    
    
    network_feature_dict=Network_Features(allnetworkspikes)
    
    
    network_feature_dict['Channel Labels']=netwrokschs
    
    network_feature_dict['Spikes per channel']=netwroksspikes
    
    network_feature_dict['Spikes All']=allnetworkspikes
    
    
    sttcs=[netpairsttcinfunc(spikeschannel, spikeschannellabel,  10, dtscale) for spikeschannel,
                                  
                                  spikeschannellabel in zip(netwroksspikes, netwrokschs)]
    
    network_feature_dict['STTC']=sttcs
    
    
    network_feature_dict['NetworkMeanSTTC']=[np.mean(np.array(l)) for l in sttcs]
    
    network_feature_dict['NetworkMedianSTTC']=[np.median(np.array(l)) for l in sttcs]
    
    network_feature_dict['Network25STTC']=[np.percentile(np.array(l), 25) for l in sttcs]
    
    
    network_df=pd.DataFrame(network_feature_dict)
    
    
    return network_df


# In[59]:


def netpairsttcinfunc(spikeschannel, spikeschannellabel,  pad, dtscale):
    
    '''takes a single bursts and returns pairwise syhnchrony between pairs
    dtscale in ms and 5ms is preferable
    pad in ms but pad does not change much but '''
    
    dict_sttc={}
    
    sttclist=[]
    
    #spikeschannel=networkburstdf['Spikes per channel']
    #spikeschannellabel=networkburstdf['Channel Labels']
    start=(sorted([sorted(i)[0] for i in spikeschannel])[0]-(pad*1000))/1000
    stop=(sorted([sorted(i)[-1] for i in spikeschannel])[-1]+(pad*1000))/1000
    
    for pair in list(combinations(range(len(spikeschannellabel)), 2)):
        
        spiketrain1 = neo.SpikeTrain(np.array(spikeschannel[pair[0]])/1000, units='ms', t_start=start,  t_stop=stop)
        spiketrain2 = neo.SpikeTrain(np.array(spikeschannel[pair[1]])/1000, units='ms', t_start=start,  t_stop=stop)
        sttc=spike_time_tiling_coefficient(spiketrain1, spiketrain2, dt=dtscale * pq.ms)
        
        dict_sttc[pair]=sttc
        sttclist.append(sttc)
        
        
    
    
    
    return sttclist


# In[8]:


def netpairsttc(networkburstdf, pad, dtscale):
    
    '''takes a single bursts and returns pairwise syhnchrony between pairs
    dtscale in ms and 5ms is preferable
    pad in ms but pad does not change much but '''
    
    dict_sttc={}
    
    sttclist=[]
    
    spikeschannel=networkburstdf['Spikes per channel']
    spikeschannellabel=networkburstdf['Channel Labels']
    start=(sorted([sorted(i)[0] for i in spikeschannel])[0]-(pad*1000))/1000
    stop=(sorted([sorted(i)[-1] for i in spikeschannel])[-1]+(pad*1000))/1000
    
    for pair in list(combinations(range(len(spikeschannellabel)), 2)):
        
        spiketrain1 = neo.SpikeTrain(np.array(spikeschannel[pair[0]])/1000, units='ms', t_start=start,  t_stop=stop)
        spiketrain2 = neo.SpikeTrain(np.array(spikeschannel[pair[1]])/1000, units='ms', t_start=start,  t_stop=stop)
        sttc=spike_time_tiling_coefficient(spiketrain1, spiketrain2, dt=dtscale * pq.ms)
        
        dict_sttc[pair]=sttc
        sttclist.append(sttc)
        
        
    networkburstdf['NetworkMeanSTTC']=np.mean(np.array(sttclist))
    
    networkburstdf['NetworkMedianSTTC']=np.median(np.array(sttclist))
    
    networkburstdf['Network25STTC']=np.percentile(np.array(sttclist), 25)
    
    
    return networkburstdf


# In[43]:


def Network_Features(networkspikes):
    
    from scipy.signal import find_peaks, peak_prominences, peak_widths
    
    
    networkdict={}
    
    
    
    networkdict['SpikeLength']=[len(l) for l in networkspikes]
    
    networkdict['Starts']=[l[0] for l in networkspikes]
    
    networkdict['Ends']=[l[-1] for l in networkspikes]
    
    networkdict['MinISI']=[min(np.diff(l)) for l in networkspikes]
    
    networkdict['MaxISI']=[max(np.diff(l)) for l in networkspikes]
    
    networkdict['MeanxISI']=[np.mean(np.diff(l)) for l in networkspikes]
    
    networkdict['MedianxISI']=[np.median(np.diff(l)) for l in networkspikes]
    
    networkdict['Rates']=[1000000/np.diff(l) for l in networkspikes]
    
    
    networkdict['Peakindexes']=[find_peaks(R)[0] for R in networkdict['Rates']]
    
    networkdict['Peakslen']=[len(peak) for peak in networkdict['Peakindexes']]
    
    networkdict['PeaksVal']=[R[indexes] for R, indexes in zip(networkdict['Rates'],  networkdict['Peakindexes'])]
    
    networkdict['PeaksValmax']=[max(pk) for pk in networkdict['PeaksVal']]
    
    
    
    networkdict['PeakProminence']=[peak_prominences(R, find_peaks(R)[0])[0] for R in networkdict['Rates']]
    
    networkdict['PeakWidth']=[peak_widths(R, find_peaks(R)[0])[0] for R in networkdict['Rates']]
    
    networkdict['MaxProm']=[max(prom) for prom in networkdict['PeakProminence']]
    
    networkdict['MaxPeakWidth']=[max(width) for width in networkdict['PeakWidth']]
    
    return networkdict


# In[15]:


def get_file_and_chid(gdf, self):
    
    ###works with meaoobj
    
    
    ###cgdf is single row for a network burst
    
    experiment=gdf['Experiment'].unique()[0]
    dlabel=gdf['Dose Label'].unique()[0]
    welllabel=gdf['Well Label'].unique()[0]
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
    
    
    
def extract_raw(gdf, self, original_fs = 20000,
        target_fs = 2000,  
        lowpass_cutoff = 500,  frequency_range1=None, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None, min_nfft=0,
                  detrend_opt='linear'):
    
    ###gdf is well and experiment and dose label groupped
    
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
    
    welldict={}
        
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
        
        burst=v[raw_start-80000:raw_start+80000].astype('float')
        
        processed_data = process_lfp(burst, original_fs, target_fs, lowpass_cutoff)
        resultsrow=process_signal(processed_data, original_fs=original_fs, target_fs=2000, chunk_size=2000, 
                                degree=10, threshold=0.10, min_duration=0.3, max_gap=0.3)
        
        welldict.update(resultsrow)
        
    return welldict
        
        
        
        
      


# In[ ]:


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
    
    
    
def extract_raw(gdf, self, original_fs = 20000,
        target_fs = 2000,  
        lowpass_cutoff = 500,  frequency_range1=None, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None, min_nfft=0,
                  detrend_opt='linear'):
    
    ###gdf is well and experiment and dose label groupped
    
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
        
        burst=v[raw_start-40000:raw_start+40000].astype('float')
        
        plt.figure()
        plt.plot(burst)
        plt.show()
        
        processed_data = process_lfp(burst, original_fs, target_fs, lowpass_cutoff)
        
        plt.figure()
        
        plt.plot(processed_data[20000:])
        plt.show()
        
        
        results=merge_regions_and_analyze_without_prominence(processed_data[20000:])
        
        welldf=pd.concat([welldf, pd.DataFrame([results])])
        
        
    return welldf


# In[ ]:


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
    
    
    
def extract_raw(gdf, self, original_fs = 20000,
        target_fs = 2000,  
        lowpass_cutoff = 500,  plotting=False):
    
    ###gdf is well and experiment and dose label groupped
    
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
        
        burst=v[raw_start-40000:raw_start+40000].astype('float')
        
        processed_data = process_lfp(burst, original_fs, target_fs, lowpass_cutoff)
        
        results=merge_regions_and_analyze_without_prominence(processed_data[20000:], plotting=plotting)
        
        welldf=pd.concat([welldf, pd.DataFrame([results])])
        
        if plotting==True:
            plt.figure()
            plt.plot(burst)
            plt.show()

            plt.figure()

            plt.plot(processed_data[20000:])
            plt.show()
        
        
    return welldf
       


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

# Function to merge regions and analyze without prominence
def merge_regions_and_analyze_without_prominence(signal, original_fs=20000, target_fs=2000, 
                                                 high_threshold=0.5, low_threshold=0.15, 
                                                 max_gap_high=0.15, max_gap_low=0.3, chunk_size=100, plotting=False):
    """
    Process signal to identify high-amplitude regions, merge regions as per thresholds and gaps, and analyze the result,
    focusing on peaks with the highest positive and negative amplitudes without using prominence.
    """
    # Identify high-amplitude regions
    
    
    from scipy.signal import find_peaks


    from scipy.signal import resample
    
    norm_signal=normalize_signal(signal, original_fs, target_fs)
    high_amplitude_regions = find_high_amplitude_regions(norm_signal, target_fs, 
                                                         threshold=high_threshold, chunk_size=chunk_size)

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
    for region in find_high_amplitude_regions(signal, fs, threshold=low_threshold, chunk_size=chunk_size):
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

    # Find maximum positive and negative peaks
    max_positive_peak_idx = np.argmax(region_signal)
    max_negative_peak_idx = np.argmin(region_signal)
    max_positive_peak = region_signal[max_positive_peak_idx]
    max_negative_peak = region_signal[max_negative_peak_idx]
    rise_rate = max_positive_peak / combined_region["Duration (seconds)"] if combined_region["Duration (seconds)"] > 0 else 0
    fall_rate = abs(max_negative_peak) / combined_region["Duration (seconds)"] if combined_region["Duration (seconds)"] > 0 else 0

    analysis_result = {
        "Start Index": combined_region["Start Index"],
        "End Index": combined_region["End Index"],
        "Duration (seconds)": combined_region["Duration (seconds)"],
        "Area Under Curve": combined_region["Area Under Curve"],
        "Max Positive Peak": max_positive_peak,
        "Max Negative Peak": max_negative_peak,
        "Rise Rate": rise_rate,
        "Fall Rate": fall_rate,
        "Positive Peak Index": max_positive_peak_idx,
        "Negative Peak Index": max_negative_peak_idx}

    # Plot the region with the selected peaks
    
    if plotting==True:
        x_values = np.arange(combined_region["Start Index"], combined_region["End Index"])
        plt.figure(figsize=(12, 6))
        plt.plot(x_values, region_signal, label="Region Signal", color="blue")
        plt.scatter(x_values[max_positive_peak_idx], region_signal[max_positive_peak_idx], color="green", label="Max Positive Peak")
        plt.scatter(x_values[max_negative_peak_idx], region_signal[max_negative_peak_idx], color="red", label="Max Negative Peak")
        plt.title("Selected Region with Extreme Peaks")
        plt.xlabel("Sample Index")
        plt.ylabel("Amplitude")
        plt.legend()
        plt.grid(True)
        plt.show()

    return analysis_result


# In[ ]:


def multitaper_spectrogram(data, fs, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None,
                           min_nfft=0, detrend_opt='linear', multiprocess=False, n_jobs=None, weighting='unity',
                           plot_on=True, return_fig=False, clim_scale=True, verbose=True, xyflip=False, ax=None):
    
    
    
    
    """ Compute multitaper spectrogram of timeseries data
    Usage:
    mt_spectrogram, stimes, sfreqs = multitaper_spectrogram(data, fs, frequency_range=None, time_bandwidth=5,
                                                            num_tapers=None, window_params=None, min_nfft=0,
                                                            detrend_opt='linear', multiprocess=False, cpus=False,
                                                            weighting='unity', plot_on=True, return_fig=False,
                                                            clim_scale=True, verbose=True, xyflip=False):
        Arguments:
                data (1d np.array): time series data -- required
                fs (float): sampling frequency in Hz  -- required
                frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] (default: [0 nyquist])
                time_bandwidth (float): time-half bandwidth product (window duration*half bandwidth of main lobe)
                                        (default: 5 Hz*s)
                num_tapers (int): number of DPSS tapers to use (default: [will be computed
                                  as floor(2*time_bandwidth - 1)])
                window_params (list): 1x2 list - [window size (seconds), step size (seconds)] (default: [5 1])
                detrend_opt (string): detrend data window ('linear' (default), 'constant', 'off')
                                      (Default: 'linear')
                min_nfft (int): minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x)
                                (default: 0)
                multiprocess (bool): Use multiprocessing to compute multitaper spectrogram (default: False)
                n_jobs (int): Number of cpus to use if multiprocess = True (default: False). Note: if default is left
                            as None and multiprocess = True, the number of cpus used for multiprocessing will be
                            all available - 1.
                weighting (str): weighting of tapers ('unity' (default), 'eigen', 'adapt');
                plot_on (bool): plot results (default: True)
                return_fig (bool): return plotted spectrogram (default: False)
                clim_scale (bool): automatically scale the colormap on the plotted spectrogram (default: True)
                verbose (bool): display spectrogram properties (default: True)
                xyflip (bool): transpose the mt_spectrogram output (default: False)
                ax (axes): a matplotlib axes to plot the spectrogram on (default: None)
        Returns:
                mt_spectrogram (TxF np array): spectral power matrix
                stimes (1xT np array): timepoints (s) in mt_spectrogram
                sfreqs (1xF np array)L frequency values (Hz) in mt_spectrogram

        Example:
        In this example we create some chirp data and run the multitaper spectrogram on it.
            import numpy as np  # import numpy
            from scipy.signal import chirp  # import chirp generation function
            # Set spectrogram params
            fs = 200  # Sampling Frequency
            frequency_range = [0, 25]  # Limit frequencies from 0 to 25 Hz
            time_bandwidth = 3  # Set time-half bandwidth
            num_tapers = 5  # Set number of tapers (optimal is time_bandwidth*2 - 1)
            window_params = [4, 1]  # Window size is 4s with step size of 1s
            min_nfft = 0  # No minimum nfft
            detrend_opt = 'constant'  # detrend each window by subtracting the average
            multiprocess = True  # use multiprocessing
            cpus = 3  # use 3 cores in multiprocessing
            weighting = 'unity'  # weight each taper at 1
            plot_on = True  # plot spectrogram
            return_fig = False  # do not return plotted spectrogram
            clim_scale = False # don't auto-scale the colormap
            verbose = True  # print extra info
            xyflip = False  # do not transpose spect output matrix

            # Generate sample chirp data
            t = np.arange(1/fs, 600, 1/fs)  # Create 10 min time array from 1/fs to 600 stepping by 1/fs
            f_start = 1  # Set chirp freq range min (Hz)
            f_end = 20  # Set chirp freq range max (Hz)
            data = chirp(t, f_start, t[-1], f_end, 'logarithmic')
            # Compute the multitaper spectrogram
            spect, stimes, sfreqs = multitaper_spectrogram(data, fs, frequency_range, time_bandwidth, num_tapers,
                                                           window_params, min_nfft, detrend_opt, multiprocess,
                                                           cpus, weighting, plot_on, return_fig, clim_scale,
                                                           verbose, xyflip):

        This code is companion to the paper:
        "Sleep Neurophysiological Dynamics Through the Lens of Multitaper Spectral Analysis"
           Michael J. Prerau, Ritchie E. Brown, Matt T. Bianchi, Jeffrey M. Ellenbogen, Patrick L. Purdon
           December 7, 2016 : 60-92
           DOI: 10.1152/physiol.00062.2015
         which should be cited for academic use of this code.

         A full tutorial on the multitaper spectrogram can be found at: # https://www.sleepEEG.org/multitaper

        Copyright 2021 Michael J. Prerau Laboratory. - https://www.sleepEEG.org
        Authors: Michael J. Prerau, Ph.D., Thomas Possidente, Mingjian He

  __________________________________________________________________________________________________________________
    """

    #  Process user input
    [data, fs, frequency_range, time_bandwidth, num_tapers,
     winsize_samples, winstep_samples, window_start,
     num_windows, nfft, detrend_opt, plot_on, verbose] = process_input(data, fs, frequency_range, time_bandwidth,
                                                                       num_tapers, window_params, min_nfft,
                                                                       detrend_opt, plot_on, verbose)

    # Set up spectrogram parameters
    [window_idxs, stimes, sfreqs, freq_inds] = process_spectrogram_params(fs, nfft, frequency_range, window_start,
                                                                          winsize_samples)
    # Display spectrogram parameters
    if verbose:
        display_spectrogram_props(fs, time_bandwidth, num_tapers, [winsize_samples, winstep_samples], frequency_range,
                                  nfft, detrend_opt)

    # Split data into segments and preallocate
    data_segments = data[window_idxs]

    # COMPUTE THE MULTITAPER SPECTROGRAM
    #     STEP 1: Compute DPSS tapers based on desired spectral properties
    #     STEP 2: Multiply the data segment by the DPSS Tapers
    #     STEP 3: Compute the spectrum for each tapered segment
    #     STEP 4: Take the mean of the tapered spectra

    # Compute DPSS tapers (STEP 1)
    dpss_tapers, dpss_eigen = dpss(winsize_samples, time_bandwidth, num_tapers, return_ratios=True)
    dpss_eigen = np.reshape(dpss_eigen, (num_tapers, 1))

    # pre-compute weights
    if weighting == 'eigen':
        wt = dpss_eigen / num_tapers
    elif weighting == 'unity':
        wt = np.ones(num_tapers) / num_tapers
        wt = np.reshape(wt, (num_tapers, 1))  # reshape as column vector
    else:
        wt = 0

    tic = timeit.default_timer()  # start timer

    # Set up calc_mts_segment() input arguments
    mts_params = (dpss_tapers, nfft, freq_inds, detrend_opt, num_tapers, dpss_eigen, weighting, wt)

    if multiprocess:  # use multiprocessing
        n_jobs = max(cpu_count() - 1, 1) if n_jobs is None else n_jobs
        mt_spectrogram = np.vstack(Parallel(n_jobs=n_jobs)(delayed(calc_mts_segment)(
            data_segments[num_window, :], *mts_params) for num_window in range(num_windows)))

    else:  # if no multiprocessing, compute normally
        mt_spectrogram = np.apply_along_axis(calc_mts_segment, 1, data_segments, *mts_params)

    # Compute one-sided PSD spectrum
    mt_spectrogram = mt_spectrogram.T
    dc_select = np.where(sfreqs == 0)[0]
    nyquist_select = np.where(sfreqs == fs/2)[0]
    select = np.setdiff1d(np.arange(0, len(sfreqs)), np.concatenate((dc_select, nyquist_select)))

    mt_spectrogram = np.vstack([mt_spectrogram[dc_select, :], 2*mt_spectrogram[select, :],
                               mt_spectrogram[nyquist_select, :]]) / fs

    # Flip if requested
    if xyflip:
        mt_spectrogram = mt_spectrogram.T

    # End timer and get elapsed compute time
    toc = timeit.default_timer()
    if verbose:
        print("\n Multitaper compute time: " + "%.2f" % (toc - tic) + " seconds")

    if np.all(mt_spectrogram.flatten() == 0):
        print("\n Data was all zeros, no output")

    # Plot multitaper spectrogram
    if plot_on:
        # convert from power to dB
        spect_data = nanpow2db(mt_spectrogram)

        # Set x and y axes
        dx = stimes[1] - stimes[0]
        dy = sfreqs[1] - sfreqs[0]
        extent = [stimes[0]-dx, stimes[-1]+dx, sfreqs[-1]+dy, sfreqs[0]-dy]

        # Plot spectrogram
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()
        im = ax.imshow(spect_data, extent=extent, aspect='auto')
        fig.colorbar(im, ax=ax, label='PSD (dB)', shrink=0.8)
        ax.set_xlabel("Time (HH:MM:SS)")
        ax.set_ylabel("Frequency (Hz)")
        im.set_cmap(plt.cm.get_cmap('Spectral_r'))
        ax.invert_yaxis()

        # Scale colormap
        if clim_scale:
            clim = np.percentile(spect_data, [5, 98])  # from 5th percentile to 98th
            im.set_clim(clim)  # actually change colorbar scale

        fig.show()
        if return_fig:
            return mt_spectrogram, stimes, sfreqs, (fig, ax)

    return mt_spectrogram, stimes, sfreqs


# Helper Functions #

# Process User Inputs #
def process_input(data, fs, frequency_range=None, time_bandwidth=5, num_tapers=None, window_params=None, min_nfft=0,
                  detrend_opt='linear', plot_on=True, verbose=True):
    """ Helper function to process multitaper_spectrogram() arguments
            Arguments:
                    data (1d np.array): time series data-- required
                    fs (float): sampling frequency in Hz  -- required
                    frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] (default: [0 nyquist])
                    time_bandwidth (float): time-half bandwidth product (window duration*half bandwidth of main lobe)
                                            (default: 5 Hz*s)
                    num_tapers (int): number of DPSS tapers to use (default: None [will be computed
                                      as floor(2*time_bandwidth - 1)])
                    window_params (list): 1x2 list - [window size (seconds), step size (seconds)] (default: [5 1])
                    min_nfft (int): minimum allowable NFFT size, adds zero padding for interpolation (closest 2^x)
                                    (default: 0)
                    detrend_opt (string): detrend data window ('linear' (default), 'constant', 'off')
                                          (Default: 'linear')
                    plot_on (True): plot results (default: True)
                    verbose (True): display spectrogram properties (default: true)
            Returns:
                    data (1d np.array): same as input
                    fs (float): same as input
                    frequency_range (list): same as input or calculated from fs if not given
                    time_bandwidth (float): same as input or default if not given
                    num_tapers (int): same as input or calculated from time_bandwidth if not given
                    winsize_samples (int): number of samples in single time window
                    winstep_samples (int): number of samples in a single window step
                    window_start (1xm np.array): array of timestamps representing the beginning time for each window
                    num_windows (int): number of windows in the data
                    nfft (int): length of signal to calculate fft on
                    detrend_opt ('string'): same as input or default if not given
                    plot_on (bool): same as input
                    verbose (bool): same as input
    """

    # Make sure data is 1 dimensional np array
    if len(data.shape) != 1:
        if (len(data.shape) == 2) & (data.shape[1] == 1):  # if it's 2d, but can be transferred to 1d, do so
            data = np.ravel(data[:, 0])
        elif (len(data.shape) == 2) & (data.shape[0] == 1):  # if it's 2d, but can be transferred to 1d, do so
            data = np.ravel(data.T[:, 0])
        else:
            raise TypeError("Input data is the incorrect dimensions. Should be a 1d array with shape (n,) where n is                             the number of data points. Instead data shape was " + str(data.shape))

    # Set frequency range if not provided
    if frequency_range is None:
        frequency_range = [0, fs / 2]

    # Set detrending method
    detrend_opt = detrend_opt.lower()
    if detrend_opt != 'linear':
        if detrend_opt in ['const', 'constant']:
            detrend_opt = 'constant'
        elif detrend_opt in ['none', 'false', 'off']:
            detrend_opt = 'off'
        else:
            raise ValueError("'" + str(detrend_opt) + "' is not a valid argument for detrend_opt. The choices " +
                             "are: 'constant', 'linear', or 'off'.")
    # Check if frequency range is valid
    if frequency_range[1] > fs / 2:
        frequency_range[1] = fs / 2
        warnings.warn('Upper frequency range greater than Nyquist, setting range to [' +
                      str(frequency_range[0]) + ', ' + str(frequency_range[1]) + ']')

    # Set number of tapers if none provided
    if num_tapers is None:
        num_tapers = math.floor(2 * time_bandwidth) - 1

    # Warn if number of tapers is suboptimal
    if num_tapers != math.floor(2 * time_bandwidth) - 1:
        warnings.warn('Number of tapers is optimal at floor(2*TW) - 1. consider using ' +
                      str(math.floor(2 * time_bandwidth) - 1))

    # If no window params provided, set to defaults
    if window_params is None:
        window_params = [5, 1]

    # Check if window size is valid, fix if not
    if window_params[0] * fs % 1 != 0:
        winsize_samples = round(window_params[0] * fs)
        warnings.warn('Window size is not divisible by sampling frequency. Adjusting window size to ' +
                      str(winsize_samples / fs) + ' seconds')
    else:
        winsize_samples = window_params[0] * fs

    # Check if window step is valid, fix if not
    if window_params[1] * fs % 1 != 0:
        winstep_samples = round(window_params[1] * fs)
        warnings.warn('Window step size is not divisible by sampling frequency. Adjusting window step size to ' +
                      str(winstep_samples / fs) + ' seconds')
    else:
        winstep_samples = window_params[1] * fs

    # Get total data length
    len_data = len(data)

    # Check if length of data is smaller than window (bad)
    if len_data < winsize_samples:
        raise ValueError("\nData length (" + str(len_data) + ") is shorter than window size (" +
                         str(winsize_samples) + "). Either increase data length or decrease window size.")

    # Find window start indices and num of windows
    window_start = np.arange(0, len_data - winsize_samples + 1, winstep_samples)
    num_windows = len(window_start)

    # Get num points in FFT
    if min_nfft == 0:  # avoid divide by zero error in np.log2(0)
        nfft = max(2 ** math.ceil(np.log2(abs(winsize_samples))), winsize_samples)
    else:
        nfft = max(max(2 ** math.ceil(np.log2(abs(winsize_samples))), winsize_samples),
                   2 ** math.ceil(np.log2(abs(min_nfft))))

    return ([data, fs, frequency_range, time_bandwidth, num_tapers,
             int(winsize_samples), int(winstep_samples), window_start, num_windows, nfft,
             detrend_opt, plot_on, verbose])


# PROCESS THE SPECTROGRAM PARAMETERS #
def process_spectrogram_params(fs, nfft, frequency_range, window_start, datawin_size):
    """ Helper function to create frequency vector and window indices
        Arguments:
             fs (float): sampling frequency in Hz  -- required
             nfft (int): length of signal to calculate fft on -- required
             frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] -- required
             window_start (1xm np array): array of timestamps representing the beginning time for each
                                          window -- required
             datawin_size (float): seconds in one window -- required
        Returns:
            window_idxs (nxm np array): indices of timestamps for each window
                                        (nxm where n=number of windows and m=datawin_size)
            stimes (1xt np array): array of times for the center of the spectral bins
            sfreqs (1xf np array): array of frequency bins for the spectrogram
            freq_inds (1d np array): boolean array of which frequencies are being analyzed in
                                      an array of frequencies from 0 to fs with steps of fs/nfft
    """

    # create frequency vector
    df = fs / nfft
    sfreqs = np.arange(0, fs, df)

    # Get frequencies for given frequency range
    freq_inds = (sfreqs >= frequency_range[0]) & (sfreqs <= frequency_range[1])
    sfreqs = sfreqs[freq_inds]

    # Compute times in the middle of each spectrum
    window_middle_samples = window_start + round(datawin_size / 2)
    stimes = window_middle_samples / fs

    # Get indexes for each window
    window_idxs = np.atleast_2d(window_start).T + np.arange(0, datawin_size, 1)
    window_idxs = window_idxs.astype(int)

    return [window_idxs, stimes, sfreqs, freq_inds]


# DISPLAY SPECTROGRAM PROPERTIES
def display_spectrogram_props(fs, time_bandwidth, num_tapers, data_window_params, frequency_range, nfft, detrend_opt):
    """ Prints spectrogram properties
        Arguments:
            fs (float): sampling frequency in Hz  -- required
            time_bandwidth (float): time-half bandwidth product (window duration*1/2*frequency_resolution) -- required
            num_tapers (int): number of DPSS tapers to use -- required
            data_window_params (list): 1x2 list - [window length(s), window step size(s)] -- required
            frequency_range (list): 1x2 list - [<min frequency>, <max frequency>] -- required
            nfft(float): number of fast fourier transform samples -- required
            detrend_opt (str): detrend data window ('linear' (default), 'constant', 'off') -- required
        Returns:
            This function does not return anything
    """

    data_window_params = np.asarray(data_window_params) / fs

    # Print spectrogram properties
    print("Multitaper Spectrogram Properties: ")
    print('     Spectral Resolution: ' + str(2 * time_bandwidth / data_window_params[0]) + 'Hz')
    print('     Window Length: ' + str(data_window_params[0]) + 's')
    print('     Window Step: ' + str(data_window_params[1]) + 's')
    print('     Time Half-Bandwidth Product: ' + str(time_bandwidth))
    print('     Number of Tapers: ' + str(num_tapers))
    print('     Frequency Range: ' + str(frequency_range[0]) + "-" + str(frequency_range[1]) + 'Hz')
    print('     NFFT: ' + str(nfft))
    print('     Detrend: ' + detrend_opt + '\n')


# NANPOW2DB
def nanpow2db(y):
    """ Power to dB conversion, setting bad values to nans
        Arguments:
            y (float or array-like): power
        Returns:
            ydB (float or np array): inputs converted to dB with 0s and negatives resulting in nans
    """

    if isinstance(y, int) or isinstance(y, float):
        if y == 0:
            return np.nan
        else:
            ydB = 10 * np.log10(y)
    else:
        if isinstance(y, list):  # if list, turn into array
            y = np.asarray(y)
        y = y.astype(float)  # make sure it's a float array so we can put nans in it
        y[y == 0] = np.nan
        ydB = 10 * np.log10(y)

    return ydB


# Helper #
def is_outlier(data):
    smad = 1.4826 * np.median(abs(data - np.median(data)))  # scaled median absolute deviation
    outlier_mask = abs(data-np.median(data)) > 3*smad  # outliers are more than 3 smads away from median
    outlier_mask = (outlier_mask | np.isnan(data) | np.isinf(data))
    return outlier_mask


# CALCULATE MULTITAPER SPECTRUM ON SINGLE SEGMENT
def calc_mts_segment(data_segment, dpss_tapers, nfft, freq_inds, detrend_opt, num_tapers, dpss_eigen, weighting, wt):
    """ Helper function to calculate the multitaper spectrum of a single segment of data
        Arguments:
            data_segment (1d np.array): One window worth of time-series data -- required
            dpss_tapers (2d np.array): Parameters for the DPSS tapers to be used.
                                       Dimensions are (num_tapers, winsize_samples) -- required
            nfft (int): length of signal to calculate fft on -- required
            freq_inds (1d np array): boolean array of which frequencies are being analyzed in
                                      an array of frequencies from 0 to fs with steps of fs/nfft
            detrend_opt (str): detrend data window ('linear' (default), 'constant', 'off')
            num_tapers (int): number of tapers being used
            dpss_eigen (np array):
            weighting (str):
            wt (int or np array):
        Returns:
            mt_spectrum (1d np.array): spectral power for single window
    """

    # If segment has all zeros, return vector of zeros
    if all(data_segment == 0):
        ret = np.empty(sum(freq_inds))
        ret.fill(0)
        return ret

    if any(np.isnan(data_segment)):
        ret = np.empty(sum(freq_inds))
        ret.fill(np.nan)
        return ret

    # Option to detrend data to remove low frequency DC component
    if detrend_opt != 'off':
        data_segment = detrend(data_segment, type=detrend_opt)

    # Multiply data by dpss tapers (STEP 2)
    tapered_data = np.multiply(np.mat(data_segment).T, np.mat(dpss_tapers.T))

    # Compute the FFT (STEP 3)
    fft_data = np.fft.fft(tapered_data, nfft, axis=0)

    # Compute the weighted mean spectral power across tapers (STEP 4)
    spower = np.power(np.imag(fft_data), 2) + np.power(np.real(fft_data), 2)
    if weighting == 'adapt':
        # adaptive weights - for colored noise spectrum (Percival & Walden p368-370)
        tpower = np.dot(np.transpose(data_segment), (data_segment/len(data_segment)))
        spower_iter = np.mean(spower[:, 0:2], 1)
        spower_iter = spower_iter[:, np.newaxis]
        a = (1 - dpss_eigen) * tpower
        for i in range(3):  # 3 iterations only
            # Calc the MSE weights
            b = np.dot(spower_iter, np.ones((1, num_tapers))) / ((np.dot(spower_iter, np.transpose(dpss_eigen))) +
                                                                 (np.ones((nfft, 1)) * np.transpose(a)))
            # Calc new spectral estimate
            wk = (b**2) * np.dot(np.ones((nfft, 1)), np.transpose(dpss_eigen))
            spower_iter = np.sum((np.transpose(wk) * np.transpose(spower)), 0) / np.sum(wk, 1)
            spower_iter = spower_iter[:, np.newaxis]

        mt_spectrum = np.squeeze(spower_iter)

    else:
        # eigenvalue or uniform weights
        mt_spectrum = np.dot(spower, wt)
        mt_spectrum = np.reshape(mt_spectrum, nfft)  # reshape to 1D

    return mt_spectrum[freq_inds]


# In[ ]:


def mtcsdfast_single_channel(x, nFFT=1024, Fs=2, WinLength=None, nOverlap=None,
                             NW=3, DetrendType='linear', FreqRange=None):
    
    
    from scipy.signal import spectrogram, get_window, detrend
    from scipy.fftpack import fft
    from scipy.signal.windows import dpss
    import matplotlib.pyplot as plt

    
    
    """
    Multitaper Spectral Density Estimation for Single Channel in Python.

    Parameters:
        x: ndarray
            Input time series, shape (nSamples,).
        nFFT: int
            Number of points for FFT computation.
        Fs: int
            Sampling frequency.
        WinLength: int
            Length of moving window.
        nOverlap: int
            Overlap between successive windows.
        NW: float
            Time-bandwidth product.
        DetrendType: str
            Type of detrending ('linear' or None).
        FreqRange: list or tuple
            Frequency range for output.

    Returns:
        psd: ndarray
            Power spectral density estimates.
        freqs: ndarray
            Frequencies at which the estimates are computed.
    """
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError("Input signal must be one-dimensional for single-channel analysis.")

    nSamples = len(x)

    if WinLength is None:
        WinLength = nFFT
    if nOverlap is None:
        nOverlap = WinLength // 2
    if FreqRange is None:
        FreqRange = (0, Fs / 2)

    # Generate Slepian sequences (tapers)
    tapers, eigenvalues = dpss(WinLength, NW, return_ratios=True)

    # Setup frequency bins
    freq_bins = np.fft.rfftfreq(nFFT, 1 / Fs)
    freq_mask = (freq_bins >= FreqRange[0]) & (freq_bins <= FreqRange[1])
    selected_freq_bins = freq_bins[freq_mask]

    # Determine number of segments
    step = WinLength - nOverlap
    num_segments = (nSamples - WinLength) // step + 1

    # Preallocate PSD
    psd = np.zeros(len(selected_freq_bins))

    for seg_idx in range(num_segments):
        start_idx = seg_idx * step
        end_idx = start_idx + WinLength
        segment = x[start_idx:end_idx]

        # Detrend if needed
        if DetrendType == 'linear':
            segment = detrend(segment, type='linear')

        # Apply tapers and compute FFT
        tapered_segments = segment[:, None] * tapers
        fft_output = fft(tapered_segments, n=nFFT, axis=0)[freq_mask]

        # Average over tapers and accumulate PSD
        psd += np.mean(np.abs(fft_output) ** 2, axis=1)

    psd /= num_segments
    
    return psd, selected_freq_bins

def plot_psd(psd, freqs):
    """Plot the Power Spectral Density."""
    plt.figure(figsize=(8, 5))
    plt.plot(freqs, psd, label="PSD", color="blue")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power Spectral Density")
    plt.title("Power Spectral Density Estimate")
    plt.grid(True)
    plt.legend()
    plt.show()




# In[17]:


def plot_spectrogram_fixed_range(data, nFFT=1024, Fs=1024, WinLength=512, 
                                 nOverlap=384, cmap='Spectral_r', FreqRange=(0, 15)):
    """
    Plot a spectrogram with a fixed frequency range and no dynamic scaling.

    Parameters:
        data: ndarray
            Input signal data.
        nFFT: int
            Number of points for FFT computation.
        Fs: int
            Sampling frequency.
        WinLength: int
            Length of the moving window.
        nOverlap: int
            Overlap between successive windows.
        cmap: str
            Colormap for the spectrogram.
        FreqRange: tuple
            Frequency range for the spectrogram display.
    """
    # Define the step size
    step = WinLength - nOverlap
    nSegments = (len(data) - WinLength) // step + 1

    # Time vector for spectrogram
    times = np.arange(nSegments) * step / Fs

    # Generate Slepian sequences (tapers)
    tapers, _ = dpss(WinLength, NW=3, return_ratios=True)

    # Prepare spectrogram matrix
    spectrogram_matrix = []

    for seg_idx in range(nSegments):
        start_idx = seg_idx * step
        end_idx = start_idx + WinLength
        segment = data[start_idx:end_idx]

        # Detrend the segment
        segment = detrend(segment, type='linear')

        # Apply tapers and compute FFT
        tapered_segment = segment[:, None] * tapers
        fft_output = fft(tapered_segment, n=nFFT, axis=0)
        power_spectrum = np.mean(np.abs(fft_output[:nFFT // 2 + 1, :])**2, axis=1)
        spectrogram_matrix.append(power_spectrum)

    spectrogram_matrix = np.array(spectrogram_matrix).T
    freqs = np.fft.rfftfreq(nFFT, d=1/Fs)

    # Limit frequency range
    freq_mask = (freqs >= FreqRange[0]) & (freqs <= FreqRange[1])
    spectrogram_matrix = spectrogram_matrix[freq_mask, :]
    freqs = freqs[freq_mask]

    # Plot the spectrogram
    plt.figure(figsize=(10, 6))
    plt.imshow(
        10 * np.log10(spectrogram_matrix),
        aspect='auto',
        origin='lower',
        extent=[times[0], times[-1], freqs[0], freqs[-1]],
        cmap=cmap
    )
    plt.colorbar(label="Power (dB)")
    plt.xlabel("Time (s)")
    plt.ylabel("Frequency (Hz)")
    plt.title("Spectrogram (Fixed Range: 0-15 Hz)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()


# Plot the spectrogram with a fixed frequency range for the new dataset
#plot_spectrogram_fixed_range(data_burstfree, nFFT=1024, Fs=1024, WinLength=512, nOverlap=384, cmap='Spectral_r', FreqRange=(0, 15))


# In[ ]:




