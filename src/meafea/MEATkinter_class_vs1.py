# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 15:23:28 2024

@author: Amalya Hakobyan, ChatGPT

to dos: 1. User channel select, by label, add ch label too. 2. User analysis segment select too.
 3, Save htsml in dpath,mkdir interacrive

"""

import pandas as pd 
import numpy as np
import h5py
import matplotlib.pyplot as plt
from sys import exit
from lxml import etree
from scipy import signal 
from sklearn import preprocessing
from scipy import optimize
import regex as re
import timeit
from lmfit.models import StepModel, LinearModel
import sys, importlib, os
from sklearn.neighbors import LocalOutlierFactor
import pickle
from sklearn.decomposition import PCA

import McsPy.McsData
import McsPy.McsCMOS
from McsPy import ureg, Q_

import statsmodels.api as sm # to build a LOWESS model
from scipy.interpolate import interp1d 
import os
import tkinter as tk
from tkinter import ttk
import webview
import threading
from bokeh.server.server import Server
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.plotting import figure, show, output_notebook
from bokeh.io import push_notebook
import nbimporter
from pathlib import Path
from module_V5Radar_vs1 import *
from moduleV15_vs1 import *
from MEAFEA_functions_vs1 import *
from Templates_vs1 import *
from BurstDetectionPossion_lib_vs1 import *

                                                                                                                                                                                                                                                                                                                                                                                                                                                            
from bokeh.io import output_notebook, show

from bokeh.models import LinearColorMapper, ColorBar, BasicTicker, FixedTicker, FuncTickFormatter
from bokeh.layouts import column, gridplot
from bokeh.models import Slider, ColumnDataSource, Span, BoxAnnotation, Select
from bokeh.plotting import figure, curdoc
from bokeh.models import ColumnDataSource, CustomJS, HoverTool, CustomJS, Div
from bokeh.layouts import gridplot
from bokeh.transform import linear_cmap
from bokeh.palettes import Viridis256

from bokeh.io import push_notebook
from ipywidgets import interact

from scipy.stats import median_abs_deviation
output_notebook()




#maybe for a single well. 


##todo 
##  vizualise spike detection and threshold per group and experiment first. 
### oad raw data, detect a threshold, show detections, plot qc. 














def choose_network_bursts(wellnetworbursts, value=-10):
    
    
    if value==-10:
        networkburstchosen=wellnetworbursts[
       (wellnetworbursts['Network_Mean_Spike_Time_Tiling_Coefficient']==wellnetworbursts['Network_Mean_Spike_Time_Tiling_Coefficient'].min()) |
       (wellnetworbursts['Network_Mean_Spike_Time_Tiling_Coefficient']==wellnetworbursts['Network_Mean_Spike_Time_Tiling_Coefficient'].max())]
    else:
        valmin=value-0.25
        valmax=value+0.25
        
        networkburstchoseninter=wellnetworbursts[
       (wellnetworbursts['Network_Mean_Spike_Time_Tiling_Coefficient']>=valmin) &
      
       (wellnetworbursts['Network_Mean_Spike_Time_Tiling_Coefficient']<=valmax)]
        
        networkburstchosen=networkburstchoseninter[
       (networkburstchoseninter['Network_Mean_Spike_Time_Tiling_Coefficient']==networkburstchoseninter['Network_Mean_Spike_Time_Tiling_Coefficient'].min()) |
       (networkburstchoseninter['Network_Mean_Spike_Time_Tiling_Coefficient']==networkburstchoseninter['Network_Mean_Spike_Time_Tiling_Coefficient'].max())]
    
    return networkburstchosen
       


def detect_single(signal,  methodarg, dead_time=0.5, threshold=5.5, negative=True, positive=True, fs=20000, mindur=3, 
                  start=0, detection=False, durationofpart=5):
    
    '''durationofpart is in second'''
    
 
    splitl=durationofpart*20000
    split=int(len(signal)//(100000)) ##10000 is half a seconf 
    parted_signal=np.array_split(signal, split)
    
    
    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*mindur
    
    


    pretrig=(fs/1000)*1
    postrig=(fs/1000)*2

    ##maximum spike duration
    
    maxdur=duration*1.5
    


    ##thresholds
    stds= [median_abs_deviation(subsignal[:]) for subsignal in parted_signal]   ##np.std##median_abs_deviation
    
    stdused=np.std(signal[:])
    
    method=methodarg[0]
    
    marg=methodarg[1]
    
    minstd=method(stds, *marg)
    
    
    meanstd=np.mean(stds)
    
    medianstd=np.median(stds)
    
    stdlist=[minstd, medianstd]
    
    if detection==True:
        
    
        if positive==False:

            negpeak_idxs, pospeak_idxs=peaknew(signal, positive, negative)



            threshold_crossings_idx=np.where(signal <(-threshold*medianstd))[0]

            #print(negpeak_idxs, pospeak_idxs, threshold_crossings_idx)

            peak_idx=np.intersect1d(negpeak_idxs, threshold_crossings_idx)

        else:
            negpeak_idxs, pospeak_idxs=peaknew(signal, positive, negative)

            threshold_crossings_idx=np.where((signal>(threshold)*medianstd) | (signal<(-threshold)*medianstd))[0]

            #print(negpeak_idxs[:10], pospeak_idxs[:10], threshold_crossings_idx[:10])

            peak_idx=negpeak_idxs+pospeak_idxs



            peak_idx=np.intersect1d(peak_idx, threshold_crossings_idx)



        ###should not start with spike
        indexes=[i for i in range(len(peak_idx)) if (((peak_idx[i]-int(duration))>0) & ((peak_idx[i]+duration)<len(signal)))]
        mins=np.array(sorted(peak_idx[indexes]))


        """removing spikes that violate dead time criteria also 50us peaks"""

        while np.any(np.diff(mins)<dead_time):

            dead_idx=(np.where(np.diff(mins)<dead_time)[0]+1).tolist()
            mins=np.delete(mins, dead_idx)




            
        
    else:
        mins=[]
        

    return ((mins), threshold*minstd)



def map_value_to_color(value, color_mapper):
    # Clip value to ensure it's within bounds
    value = np.clip(value, color_mapper.low, color_mapper.high)
    
    # Normalize the value to a 0-1 range
    norm_value = (value - color_mapper.low) / (color_mapper.high - color_mapper.low)
    
    # Map the normalized value to an index in the color palette
    color_index = int(norm_value * (len(color_mapper.palette) - 1))
    
    # Return the color from the palette
    return color_mapper.palette[color_index]

def time_reverse(x, start):
    
    x=(x-start)/50
    
    return x

def choose_channel_for_plot(welldf):
    
    
    dictchannels={}
    
    dictchannels['MinNormal']=welldf[welldf['lof_scores']==welldf['lof_scores'].min()]['Channel ID'].unique()[0]
    
    dictchannels['MaxNormal']=welldf[welldf['lof_scores']==welldf['lof_scores'].max()]['Channel ID'].unique()[0]
    
    welldfmed=welldf[welldf['lof_scores']>=np.median(welldf['lof_scores'].values)].sort_values('lof_scores')
    dictchannels['MedianNormal']=welldfmed['Channel ID'].unique()[0]
    
    
    return  dictchannels 



def outlier_channel(combinedbychanneldata):
    
    
    from sklearn.neighbors import LocalOutlierFactor
    
    bdatf=combinedbychanneldata[['MFR(Hz)', 'BFR(Hz)', 'Q1 Burst_Duration_(sec)', 'CV Number_of_Spikes_per_Burst', 
                     'Mean Burstiness_Index', 'Mean Burst_Variance_of_Spike_Times','Mean IBI','Q1 IBI','CV ISI']]
    
    lof = LocalOutlierFactor(n_neighbors = 20)
    #fit it to the training data, since we don't use it for novelty than this is fine
    y_pred = lof.fit_predict(bdatf)
    #extract the predictions as strings
    combinedbychanneldata["lof_outliers"] = y_pred.astype(str)
    #print the number of outliers relative to non-outliers
    print(combinedbychanneldata["lof_outliers"].value_counts())
    #extract the outlier scores
    combinedbychanneldata["lof_scores"] = lof.negative_outlier_factor_
    
    return combinedbychanneldata


def outlier_burst(burstdata):
    
    bdatf=burstdata[['Burst_Duration_(sec)',
       
       'Burst_Variance_of_Spike_Times', 'Burst_Max_Inter-Spike_Interval','Burst_Mean_Inter-Spike_Interval', 
                     'Burstiness_Index']]
    
    from sklearn.ensemble import IsolationForest
    #create the method instance
    isf = IsolationForest(n_estimators = 100, random_state = 42, contamination = 0.02)
    #use fit_predict on the data as we are using all the data
    preds = isf.fit_predict(bdatf)
    #extract outliers from the data
    burstdata["iso_forest_outliers"] = preds
    burstdata["iso_forest_outliers"] = burstdata["iso_forest_outliers"].astype(str)
    #extract the scores from the data in terms of strength of outlier
    burstdata["iso_forest_scores"] = isf.decision_function(bdatf)
    
    #outlierdat=burstdata[burstdata["iso_forest_outliers"]==(-1)]
    
    
    return burstdata

def get_processed_bursts(burstdata):
    
    '''get bursts from a channel label s according to normality'''
    
    
    
    
    burstdata=outlier_burst(burstdata)
    
    color_mapper = LinearColorMapper(palette=Viridis256, 
                                     low=burstdata["iso_forest_scores"].min(), 
                                     high=burstdata["iso_forest_scores"].max())
    #color = map_value_to_color(value, color_mapper)
    
    return burstdata, color_mapper





def choose_bursts(burstdatawithoutliers, channel):
    
   
    
    burstsofchannels=burstdatawithoutliers[burstdatawithoutliers['Channel ID']==channel].sort_values(by="iso_forest_scores")
    
    medianbursts=burstsofchannels[burstsofchannels["iso_forest_scores"]>=np.median(burstsofchannels["iso_forest_scores"].values)].iloc[:5, :]
    
    maxbursts=burstsofchannels[burstsofchannels["iso_forest_scores"]<=burstsofchannels["iso_forest_scores"].max()].iloc[-5:-1, :]
    
    
    minbursts=burstsofchannels[burstsofchannels["iso_forest_scores"]>=burstsofchannels["iso_forest_scores"].min()].iloc[:5, :]
    
    return pd.concat([minbursts, medianbursts, maxbursts ], axis=0)


def PoissonSurpriseGraph(sp, f, surpriseth, failfactor, startfactor):
    
    
    ##print(surpriseth, 'supriseth')
    
    """Function to  detect burst-like spike trains"""
    
    Bursts=[] #container
   
   
    sp=np.sort(sp) #sorted spikes
    
    ISI=np.diff(sp) #ISIs 
    ISImean=np.mean(ISI)
    ISIhalfmean=ISImean/2
    ISIdoublemean=2*ISImean
    
    test_surprise=[]
    
    test_Nofspikes=[]
    
    test_spikingrate=[]
    
    #figevent, axesevent=plt.subplots(1, 1, figsize=(30, 6))
    #axesevent.eventplot(sp, colors='grey', linewidths=0.8, alpha=0.6, linelengths=0.8)
   
    
    
    
    
    
    
    #f=len(sp)/((sp[-1]-sp[0])/1000000) #N/seconds
    
    #esimench=dynsurprise(sp, ISI, f)
    
    ###print(f, 'Average FIring rate')
    start_indexes=np.where(ISI<ISIhalfmean)[0] ##strict threshold
    
    end_index=-1
    
    for i in range(len(start_indexes)-3):
        
        if  start_indexes[i]>end_index:
            bind=[-1, -1]

            if all(np.diff(start_indexes[i:i+3])==1):

                bind[0]=start_indexes[i]
                bind[1]=start_indexes[i+2] #or 4
                

                p0=Surprise(sp, bind, ISI, f)
                
                
                fwd_failed=0
                
                test_bind=bind

                for j in range(start_indexes[i+2]+1, len(ISI)):
                    
                   

                    if ISI[j]<ISIdoublemean and fwd_failed<failfactor:
                        test_bind[1]=j
                        

                        
                        p=Surprise(sp, test_bind, ISI, f)
                        
                       
                        
                        #if extension maximizes the surprise, extend

                        if p>=p0:
                            bind[1]=j
                            p0=p
                            
                        #if extension is not maximizing, bursts right index is defined
                        
                            
                       
                        else:
                            fwd_failed=fwd_failed+1
                            end_index=bind[1]
                            
                            
                          
                            #search for left index
                    else:
                        end_index=bind[1]
                        
                        
                        for i in range(bind[0]+1, bind[1]):
                            
                          
                            
                            p=Surprise(sp, [i, bind[1]], ISI, f)
                            #if removing from the start only increases surprise, remove 
                            if p>startfactor*p0:
                                bind[0]=i
                                p0=p
                            #else left is defined
                            else:
                                #if defined burst is higher than the given threshold
                                if p0>surpriseth:
                                    
                                    
                                    
                                    test_surprise.append(p0)
                                    test_Nofspikes.append(len(sp[bind[0]:bind[1]+1]))
                                    test_spikingrate.append(len(sp[bind[0]:bind[1]+1])/((sp[bind[0]:bind[1]+1][-1])-
                                                            (sp[bind[0]:bind[1]+1][0])/1000000))
                                    ###print(p, p0, 'p, p0')
                                    
                                   

                                    #axesevent.eventplot(sp[bind[0]:bind[1]+1].tolist(), colors=color, linewidths=0.8, linelengths=0.8)
                                    #plt.annotate(str(p0), xy=(sp[bind[0]:bind[1]+1].tolist()[0], 1.5), xytext=(sp[bind[0]:bind[1]+1].tolist()[0], 1.5))
                                    
                                    
                                   
                                        
                                        
                                    Bursts.append(sp[bind[0]:bind[1]+1].tolist())

                                    bind=[-1, -1]
                                else:
                                    bind=[-1, -1]

                                break
                        break 
    
    
    
  
    
    
    #plot_burst(test_surprise, test_Nofspikes, test_spikingrate)

           

    return (Bursts, test_surprise)
            
            

def PoissonSurpriseGraphMerge(sp, f, surpriseth, failfactor, startfactor, minIBI):
    
    
    ##print(surpriseth, 'supriseth')
    
    """Function to  detect burst-like spike trains"""
    
    Bursts=[] #container
    merge=[]
   
   
    sp=np.sort(sp) #sorted spikes
    plots=0
    ISI=np.diff(sp) #ISIs 
    ISImean=np.mean(ISI)
    ISIhalfmean=ISImean/2
    ISIdoublemean=2*ISImean
    
    test_surprise=[]
    
    test_Nofspikes=[]
    
    test_spikingrate=[]
   
    start_indexes=np.where(ISI<ISIhalfmean)[0] ##strict threshold
    
    end_index=-1
    
    for i in range(len(start_indexes)-3):
        
        if  start_indexes[i]>end_index:
            bind=[-1, -1]

            if all(np.diff(start_indexes[i:i+3])==1):

                bind[0]=start_indexes[i]
                bind[1]=start_indexes[i+2] #or 4
                

                p0=Surprise(sp, bind, ISI, f)
                
                
                fwd_failed=0
                
                test_bind=bind

                for j in range(start_indexes[i+2]+1, len(ISI)):
                    
                   

                    if ISI[j]<ISIdoublemean and fwd_failed<failfactor:
                        test_bind[1]=j
                        

                        
                        p=Surprise(sp, test_bind, ISI, f)
                        
                       
                        
                        #if extension maximizes the surprise, extend

                        if p>=p0:
                            bind[1]=j
                            p0=p
                            
                        #if extension is not maximizing, bursts right index is defined
                        
                            
                       
                        else:
                            fwd_failed=fwd_failed+1
                            end_index=bind[1]
                            
                            
                          
                            #search for left index
                    else:
                        end_index=bind[1]
                        
                        
                        for i in range(bind[0]+1, bind[1]):
                            
                          
                            
                            p=Surprise(sp, [i, bind[1]], ISI, f)
                            #if removing from the start only increases surprise, remove 
                            if p>startfactor*p0:
                                bind[0]=i
                                p0=p
                            #else left is defined
                            else:
                                #if defined burst is higher than the given threshold
                                if p0>surpriseth:
                                    
                                    
                                    
                                    test_surprise.append(p0)
                                    test_Nofspikes.append(len(sp[bind[0]:bind[1]+1]))
                                    test_spikingrate.append(len(sp[bind[0]:bind[1]+1])/((sp[bind[0]:bind[1]+1][-1])-
                                                            (sp[bind[0]:bind[1]+1][0])/1000000))
                                    ###print(p, p0, 'p, p0')
                                    
                                   

                                    #axesevent.eventplot(sp[bind[0]:bind[1]+1].tolist(), colors=color, linewidths=0.8, linelengths=0.8)
                                    #plt.annotate(str(p0), xy=(sp[bind[0]:bind[1]+1].tolist()[0], 1.5), xytext=(sp[bind[0]:bind[1]+1].tolist()[0], 1.5))
                                    
                                    
                                   
                                        
                                        
                                    Bursts.append(sp[bind[0]:bind[1]+1].tolist())

                                    bind=[-1, -1]
                                else:
                                    bind=[-1, -1]

                                break
                        break 
    
    
    
  
    
    
      ##print(len(Bursts))        
    if len(Bursts) > 0:       
        merge.append(Bursts[0])
        
      

        for B in Bursts[1:]:
            
            if B[0]-merge[-1][-1] < minIBI: #merge burst that are not separated with minIBI
                merge[-1].extend(B)
                
               

            else:
                merge.append(B)

           

    return (merge, test_surprise)


def MaxInterval(sp, f, startISI, endISI, minIBI, minNSpikes, minDur):
    
    """Function to  detect burst-like spike trains
    constructed accoring https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4969396/"""
    
   
    
    Bursts=[]
    
    bind=[]
    merge=[]
    
    ISI=np.diff(sp)
    ##print(min(ISI), 'minISI')
    
    start_indexes=np.where(ISI<startISI)[0]
    ##print(len(start_indexes), 'lenstarts')
    end_index=-1
    
    for start_index in start_indexes:
        
        bind=[]
        
        ###print(start_index, 'start_index')
        
        if start_index>end_index:
            
            ###print('here')
        
            for i in range(start_index, len(ISI)):
                
                ###print(i, ISI[i], 'isi')
                if ISI[i]<endISI:
                    bind.append(i)
                else:
                    end_index=i
                    
                    if (len(bind)>=minNSpikes) and (sp[max(bind)+1]-sp[min(bind)])>=minDur:
                        
                        Bursts.append(sp[min(bind):(max(bind)+1)].tolist())
                        bind=[]
                        
                    bind=[]
                    break
                    
        else:
            continue
            
     
                    
    ##print(len(Bursts))        
    if len(Bursts) > 0:       
        merge.append(Bursts[0])
        
      

        for B in Bursts[1:]:
            
            if B[0]-merge[-1][-1] < minIBI: #merge burst that are not separated with minIBI
                merge[-1].extend(B)
                
               

            else:
                merge.append(B)
                
                
    return (merge, np.ones(len(merge)))
            
            



def burst_and_features(spikes, duration, function, args):
    
    df=pd.DataFrame()
    
    spiketrain=spikes['Timestamp [µs]'].sort_values(ascending=True).values
    f=len(spiketrain)/duration
    bursts=function(spiketrain, f, *args)[0]
    
    burstspikes=sum([len(bur) for bur in bursts])
        
    ISI=np.diff(spiketrain)
    
    for burst in bursts:
        
        if len(burst)>0:
            
            s=Surprise(spiketrain, burst, ISI, f, external=True)
            b=np.array(burst)
            tempdf=pd.DataFrame(spikes.iloc[0, 0:8]).T
           
            tempdf['Start timestamp [µs]']=b[0]
            tempdf['Duration']=(np.log10(b[-1]-b[0]))
            tempdf['Spike Count']=len(b)

            tempdf['FT Spikes']=len(b)*100/len(spiketrain)

            tempdf['Max ISI']=max(np.log10(np.diff(b)))
            tempdf['Min ISI']=min(np.log10(np.diff(b)))
            tempdf['Mean ISI']=np.mean(np.log10(np.diff(b)))

            tempdf['Variance']=np.std(np.log10(np.diff(b)))

            tempdf['Surprise']=s
            tempdf['Burstiness']=burstspikes/len(spiketrain)



            #tempdf['Spikes']=tempdf['Spikes'].astype('object')

            tempdf.loc[:, 'Spikes']=[b]

           

            df=pd.concat([df, tempdf], axis=0)
            
            
    return df



def update_and_save_layout(filename, plot_list, layoutName):
    
    
    # Combine all plots into a single layout
    layout = column(plot_list)

    # Generate HTML for the combined layout
    html = file_html(layout, CDN, layoutName)

    # Write the HTML to the file
    with open(filename, 'w') as file:
        file.write(html)
        
    curdoc().clear()

def show_plot(filename, layoutName):
        # Create a web view to display the Bokeh plots
        webview.create_window(layoutName, filename)
        webview.start()
    
          



class MEATkinterApp:
    
    
    def __init__(self, meaobject, comparevar):
        
        
        self.hdf5path=meaobject.hdf5path
        
        self.stamp=meaobject.stamp
        
        self.spikeparams_set=meaobject.spikeparams_set
        self.burstparams_set=meaobject.burstparams_set
        
        meaobject.spikedata['Dose Label']=meaobject.spikedata['Dose Label'].values.astype('str')
        
        self.spikedata=meaobject.spikedata
        
        self.burstdata=meaobject.burstdata
        
        self.sourcerawinteractive={}
        
        
        
        
        burstwithoutliers, color_mapper=get_processed_bursts(self.burstdata)
        self.burstwithoutliers=burstwithoutliers
        self.burstwithoutlierscolor_mapper=color_mapper
        
        self.networkburstdata=meaobject.networkburstdata
        
        self.combinedbychanneldata=meaobject.CombinedbyChannel
        #print(self.combinedbychanneldata.shape)

        self.wells=self.spikedata['Well Label'].unique()
        
        self.experiments=self.spikedata['Experiment'].unique()
        
        self.function_vars=['Default', 'Detected',
        'MaxInterval',
        'PoissonSurpriseGraph', 'PoissonSurpriseGraphMerge']
        
        self.dose_vars=self.spikedata[comparevar].unique()
        
        self.comparevar=comparevar
        
       
        
    def tkstart(self):
        
        root=tk.Tk()
    
        root.title("Data Visualization Platform")
        
       
        startisi, endisi, ibi, ns, lenb=mymea.burstparams_set['burst detection kwargs']['maxinterval']['set7']
        thresh=self.spikeparams_set['spike detection kwargs']['threshold']
        
       
        thresholdlabel = tk.Label(root, text="Enter a threshold: ")
        thresholdlabel.pack()
    
        # Entry widget to take user input
        
        
        self.thresholdentry = tk.Entry(root, textvariable=tk.StringVar(value=str(thresh)),  fg="black", bg="white")
        self.thresholdentry.pack()
        
        self.arg1label = tk.Label(root, text="Enter a Start ISI:")
        self.arg1label.pack()
    
        # Entry widget to take user input
        
        
        self.arg1entry = tk.Entry(root, textvariable=tk.StringVar(value=str(startisi)),  fg="black", bg="white" )
        self.arg1entry.pack()
    
        self.arg2label = tk.Label(root, text="Enter an END ISI:")
        self.arg2label.pack()
    
        # Entry widget to take user input
        self.arg2entry = tk.Entry(root, textvariable=tk.StringVar(value=str(endisi)),  fg="black", bg="white")
        self.arg2entry.pack()
    
        self.arg3label = tk.Label(root, text="Enter a min interburst duration:")
        self.arg3label.pack()
    
        # Entry widget to take user input
        self.arg3entry = tk.Entry(root, textvariable=tk.StringVar(value=str(ibi)),  fg="black", bg="white")
        self.arg3entry.pack()
    
        self.arg4label = tk.Label(root, text="Enter a min number of spikes:")
        self.arg4label.pack()
    
        # Entry widget to take user input
        self.arg4entry = tk.Entry(root, textvariable=tk.StringVar(value=str(ns)),  fg="black", bg="white")
        self.arg4entry.pack()
        
        
        self.arg6label = tk.Label(root, text="Enter a min duration of bursts:")
        self.arg6label.pack()
    
        # Entry widget to take user input
        self.arg6entry = tk.Entry(root, textvariable=tk.StringVar(value=str(lenb)),  fg="black", bg="white")
        self.arg6entry.pack()
    
        self.arg5label = tk.Label(root, text="Enter a Surprise: ")
        self.arg5label.pack()
    
        # Entry widget to take user input
        self.arg5entry = tk.Entry(root, textvariable=tk.StringVar(value=str(7)))
        self.arg5entry.pack()
        
        self.arg7label = tk.Label(root, text="Enter the STTC value for networks: ")
        self.arg7label.pack()
    
        # Entry widget to take user input
        
        sttc=tk.StringVar(value='-10')
        self.arg7entry = tk.Entry(root, textvariable=sttc)
        self.arg7entry.pack()
        
        self.tracedurationentry = tk.Label(root, text="Enter the Trace duration: ")
        self.tracedurationentry.pack()
    
        
        tracedur=tk.StringVar(value='10')
        self.tracedurationentry=tk.Entry(root, textvariable=tracedur)
        self.tracedurationentry.pack()
    
        # Create a slider for the X-axis range
        slider_x_label = ttk.Label(root, text="End Failures")
        slider_x_label.pack()
        self.slider_x = tk.Scale(root, from_=1, to=15, orient=tk.HORIZONTAL)
        self.slider_x.pack()
    
        # Create a slider for the Y-axis range
        slider_y_label = ttk.Label(root, text="Start Failures")
        slider_y_label.pack()
        self.slider_y = tk.Scale(root, from_=1, to=10, orient=tk.HORIZONTAL)
        self.slider_y.pack()
        
        
        well_label = ttk.Label(root, text='Select well')
        well_label.pack()
    
        self.selected_well = tk.StringVar()
        self.selected_well.set(self.wells[-1])
        
        
    
        option_menu = tk.OptionMenu(root, self.selected_well, *self.wells)
        option_menu.pack()
        
        experiment_label = ttk.Label(root, text="Select Experiment")
        experiment_label.pack()
    
        self.selected_experiment = tk.StringVar()
        
        self.well_experiments=self.spikedata[self.spikedata['Well Label']==self.selected_well.get()]['Experiment'].unique()
        self.selected_experiment.set(self.well_experiments[0])
    
        expoption_menu = tk.OptionMenu(root, self.selected_experiment, *self.well_experiments)
        expoption_menu.pack()
        
        self.ref_experiment=tk.StringVar()
        self.ref_experiment.set(self.well_experiments[0])
        
        refexpoption_menu = tk.OptionMenu(root, self.ref_experiment,
                                          *self.well_experiments)
        refexpoption_menu.pack()
        
        
        
    
        # Create a dropdown menu for the plot type
        function_label = ttk.Label(root, text="Burst Detection")
        function_label.pack()
    
        self.function_var = tk.StringVar()
        self.function_var.set('Default')
        function_menu = ttk.OptionMenu(root, self.function_var, *self.function_vars)
        function_menu.pack()
    
        dose_label = ttk.Label(root, text="Dose")
        dose_label.pack()
    
        self.dose_var = tk.StringVar()
        
        print(self.selected_experiment.get(), 'experiment')
        
        self.well_doses=self.spikedata[(self.spikedata['Well Label']==self.selected_well.get()) &
                                       (self.spikedata['Experiment']==self.selected_experiment.get())][self.comparevar].unique()
        
        
        self.dose_var.set(self.well_doses[0])
        dose_menu = ttk.OptionMenu(root, self.dose_var, *self.well_doses)
        dose_menu.pack()
        
        
        
        
        
        filter_button = ttk.Button(root, text="Apply Filters", command=self.update_filtered_data)
        filter_button.pack()
    
        #load_raw_button = tk.Button(root, text="Load raw Dataset", command=lambda : app.graph_raw(meaobject, thresholdentry.get(), filename))
        #load_raw_button.pack(pady=20)
    
        # Process Data Button
        load_burst_button = tk.Button(root, text="Load Burst Data", 
                                      command=lambda : self.graph_bursts())
        load_burst_button.pack()
        
        # Process Data Button
        load_netburst_button = tk.Button(root, text="Load Network Burst", 
                                         command=lambda : self.graph_networks())
        load_netburst_button.pack()
    
        #load_networkburst_button = tk.Button(root, text="Load Network Burst Data", command=lambda : graph_networkbursts())
        #load_networkburst_button.pack(pady=20)
        
        quit_button = tk.Button(root, text="Quit", command=root.quit)
        quit_button.pack(pady=20)
        
        root.mainloop()
        
        
        


    def update_filtered_data(self):
        
        ##generate a file nname for a bokeh plot
        
        self.filenamebursts='Bursts'+self.selected_experiment.get()+self.selected_well.get()+self.dose_var.get()+self.function_var.get()+'.html'
        
        self.filenamenets='NetworkBursts'+self.selected_experiment.get()+self.selected_well.get()+self.dose_var.get()+self.function_var.get()+'.html'
        
        ##generate a filter dict for a well level dataset selection 
        
        condition_dict={'Well Label':self.selected_well.get(), 
                        'Experiment':self.selected_experiment.get(),
                        self.comparevar:self.dose_var.get()}
        ref_condition_dict={'Well Label':self.selected_well.get(), 
                        'Experiment':self.ref_experiment.get(),
                        self.comparevar:self.dose_var.get()}
        
        
        self.wellspikes=condition_choose_query(self.spikedata, condition_dict)
        
        self.dfoutliers=condition_choose_query(self.burstwithoutliers, 
                                               condition_dict).reset_index()
        
        self.wellnetworkburstdata=condition_choose_query(self.networkburstdata, condition_dict).reset_index()
        
        
        ##channel outliers are in a experiment level
        
        
        combinedbychanneldat=self.combinedbychanneldata[self.combinedbychanneldata['Experiment']==self.selected_experiment.get()].dropna()
        combinedbychannelref=self.combinedbychanneldata[self.combinedbychanneldata['Experiment']==self.ref_experiment.get()].dropna()
        cmbychanneloutlier=outlier_channel(combinedbychanneldat)
        refcmbychanneloutlier=outlier_channel(combinedbychannelref)
        
        
        
        self.rawfile=self.wellspikes['File'].unique()[0][:-10]
        
        self.welldf=condition_choose_query(cmbychanneloutlier, condition_dict)
        self.refwelldf=condition_choose_query(refcmbychanneloutlier, 
                                              ref_condition_dict)
        
        self.dictofchannels=choose_channel_for_plot(self.refwelldf)
        
        self.mediansttc=float(self.arg7entry.get())
        self.traceduration=1000000*float(self.tracedurationentry.get())
        
    
        self.burstplot_list = []
        
        self.netplot_list=[]
        
        self.rawplot_list=[]
         
        self.plot_dict={}
        
        sourcerawinteractive={}
        
        if self.function_var.get()=='PoissonSurpriseGraph':
            
            sourcerawinteractive['args']=[int(self.arg5entry.get()), int(self.slider_x.get()), int(self.slider_y.get())]
            
            print(sourcerawinteractive['args'], 'poisson args')
            
            sourcerawinteractive['function']=PoissonSurpriseGraph
            
            
        elif self.function_var.get()=='MaxInterval': 
            
            sourcerawinteractive['args']=[int(self.arg1entry.get()), int(self.arg2entry.get()), int(self.arg3entry.get()), 
                                          int(self.arg4entry.get()), int(self.arg6entry.get())]
            
            sourcerawinteractive['function']=MaxInterval
            
        elif self.function_var.get()=='PoissonSurpriseGraphMerge':
            
            sourcerawinteractive['args']=[int(self.arg5entry.get()), int(self.slider_x.get()), int(self.slider_y.get()),
                                          int(self.arg3entry.get())]
            
            print(sourcerawinteractive['args'], 'poisson args')
            
            sourcerawinteractive['function']=PoissonSurpriseGraphMerge
            
        else:
            sourcerawinteractive={}
            
        self.sourcerawinteractive=sourcerawinteractive
        
        
    def graph_bursts(self):
        
        '''comes after update sources, runs plot channels '''
        
        if len(self.sourcerawinteractive)==0:
            
            ##this is the case when burst detection comes from the file
            
            
            ###dict of channels is based outlierness for
        
            
            
            for channelid in self.dictofchannels.values():
                
                spikes, theshold, fv, dose_st, dose_end=self.tk_get_raw(channelid)
                
                channelbursts=self.dfoutliers[self.dfoutliers['Channel ID']==channelid]
                
                #self.tk_get_raw_traces_withBursts( channelid, fv, channelbursts,
                                             #dose_st, theshold, iso=True)
                
                
                self.tk_get_raw_traces_withBurstssegment(channelid, fv, channelbursts,
                                             dose_st, dose_end, theshold, iso=True)
                
                
        else:
            print('updating the plot')
            
            function=self.sourcerawinteractive['function']
            
            args=self.sourcerawinteractive['args']
            
            for channelid in self.dictofchannels.values():
                
                spikes, theshold, fv, dose_st, dose_end=self.tk_get_raw(channelid)
                
                duration=(dose_end-dose_st)/1000000
                
                spikedf=self.wellspikes[self.wellspikes['Channel ID']==channelid]
                
                print(spikedf.shape, 'spikedf', duration, 'duration')
                
                channelbursts=burst_and_features(spikedf, duration, function, 
                                                 args).reset_index(drop=True)
                
                if len(channelbursts)>0:
                
                    channelbursts['Timestamp [µs]']=channelbursts['Start timestamp [µs]'].values
                
                    self.channelbursts=channelbursts
                
                    print(channelbursts.shape, 'channelbursts shape')
                    ###self.tk_get_raw_traces_withBursts(channelid, fv, channelbursts,
                                             #dose_st, theshold, iso=False)
                
                
                    self.tk_get_raw_traces_withBurstssegment(channelid, fv, 
                                                channelbursts,
                                                dose_st, dose_end, theshold, iso=False)
                    
                else:
                    print('No Bursts Detected on '+str(channelid))
                
                
        update_and_save_layout(self.filenamebursts,  self.burstplot_list, 'Bursts')
        
        #### Open the updated plot in a web view
        ###threading.Thread(target=show_plot, args=(self.filename, 'Bursts')).start()
        
        
    def graph_networks(self):
        
        #choose_network bursts
        
        nets=choose_network_bursts(self.wellnetworkburstdata, value=self.mediansttc)
        
        print(nets.shape, 'network shape')
        
        ##for each network on plot? as channels may be different
        
       
        
        for index, net in nets.iterrows():
            
            netlist=[]
            
            labelschannels=net['Channel Labels'][:3]
            
            channelsspikes=net['Spikes per channel'][:3]
            
            channelids=[self.wellspikes[
            self.wellspikes['True Label']==chlabel]['Channel ID'].unique()[0] 
                for chlabel in labelschannels]
            
            spk=net['Spikes All']
            
            spikst, theshold, fv, dose_st, dose_end=self.tk_get_raw(channelids[0])
            
            tms=sorted([int(time_reverse(x, dose_st)) for x in spk])
            
            
            rawst=tms[0]-40000
            rawend=tms[-1]+40000
            
            for spikes, chid, chlab in zip(channelsspikes, channelids, labelschannels):
                
                pnet=figure(width=1000, height=200, 
                            title=str(chid)+' '+str(chlab)+' '
                            +str(net['Network_Mean_Spike_Time_Tiling_Coefficient']))
                
                spikst, theshold, fv, dose_st, dose_end=self.tk_get_raw(chid)
                
                color='red'

                times=[int(time_reverse(x, dose_st)) for x in spikes]
            
                boxleft=times[0]
                boxright=times[-1]
                
                rawst=times[0]-40000
                rawend=times[-1]+40000
                
                rawsxs=np.arange(rawst, rawend, 1)
        
                rawsys=fv[int(rawst):int(rawend)]
        
        

                pnet.line(rawsxs, rawsys, color='grey')
                
                hline1 = Span(location=theshold, dimension='width', line_color='red', line_width=2, line_dash='dotted')
                        
                hline2 = Span(location=theshold*(-1), dimension='width', line_color='red', line_width=2, line_dash='dotted')
                
                pnet.add_layout(hline1)
        
                pnet.add_layout(hline2)
        
                pnet.yaxis.axis_label ='ADC'
                
                box=BoxAnnotation(left=boxleft, right=boxright, fill_alpha=0.5, 
                                  fill_color=color)
            
                pnet.add_layout(box)
                
                netlist.append(pnet)
                
            self.netplot_list.extend(netlist)  
        
        update_and_save_layout(self.filenamenets,  self.netplot_list, 'NetworkBursts')
        
        
    def tk_Raw_Threshold(self, channelid, methodarg, methodstring, argumentsnew={}, durationofpart=5):
        
        '''Single channel raw data source '''
    
        
        burstdict=choose_bursts(self.dfoutliers, channelid)
        
        channelbursts=self.dfoutliers[self.dfoutliers['Channel ID']==channelid]
        
        spikes, theshold, fv, dose_st, dose_end=self.tk_get_raw(self.file, channelid, self.dose_var.get())
        
        self.tk_get_raw_traces(channelid, fv, burstdict, dose_st,  theshold, methodstring)
        
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
        
        
    def tk_get_raw_traces_withBursts(self, channelid, fv, channelbursts, dose_st,  threshold,  iso=False):
        
        ##self, channelid, fv, bursts, dose_st,  threshold, methodstring
    
        rawsxs=[]
        rawsys=[]
        boxesleft=[]
        boxesright=[]
        shadowcolor=[]
        
        ##choosing min, median, max bursts 
        
        burstdict=choose_bursts(channelbursts, channelid)
       
            
        p=figure(width=500, height=200, title=str(channelid))
        
        appender=0
        
        
        for inindex, inrow in burstdict.iterrows():
            
            spkin=inrow['Spikes']
            
            
            ##search for the burst within -2 +2 second window
            
            searchinrangestart=spkin[0]-(50*40000)

            searchinrangeend=spkin[-1]+(50*40000)
            ##if any spike burst is in the range of this window include that bursts 
            
            channelburstsin=channelbursts[channelbursts['Spikes'].apply(lambda x: len(
                [c for c in x if ((c>=searchinrangestart) & (c<=searchinrangeend))==True])>0)]
            
            for index, row in channelburstsin.iterrows():
                
                spk=row['Spikes']
                
                if iso==True:

                    color=map_value_to_color(row['iso_forest_scores'], 
                                             self.burstwithoutlierscolor_mapper)
                
                else:
                    color='blue'

                times=[int(time_reverse(x, dose_st)) for x in spk]

                start=times[0]-40000

                if start<=0:

                    start=0


                end=times[-1]+40000



                xval=np.arange(start, end, 1)
                
                xval=xval-start+appender
                
                boxleft=times[0]-start+appender+10000
                boxright=times[-1]-start+appender+10000
                
                appender=xval[-1]
                
                xval=xval+10000

                rawsxs.extend(xval.tolist())

                yval=fv[int(start):int(end)]

                rawsys.extend(yval.tolist())
                
                boxesleft.append(boxleft)
                boxesright.append(boxright)
                shadowcolor.append(color)

        p.line(rawsxs, rawsys, color='grey')
        
        hline1 = Span(location=threshold, dimension='width', line_color='red', line_width=2, line_dash='dotted')
                
        hline2 = Span(location=threshold*(-1), dimension='width', line_color='red', line_width=2, line_dash='dotted')
        
        p.add_layout(hline1)
        
        p.add_layout(hline2)
        
        p.yaxis.axis_label ='ADC'
        
        for left, right, color in zip(boxesleft, boxesright,  shadowcolor):
            box=BoxAnnotation(left=left, right=right, fill_alpha=0.5, fill_color=color)
            
            p.add_layout(box)
        
        ###each channel is a single plot for a well.####
        
        self.burstplot_list.append(p)
        
        
    def tk_get_raw_traces_withBurstssegment(self, channelid, fv,
                                            channelbursts, dose_st, dose_end, 
                                            threshold,  iso=False):
        
       
        ##self, channelid, fv, bursts, dose_st,  threshold, methodstring
    
        rawsxs=[]
        rawsys=[]
        boxesleft=[]
        boxesright=[]
        shadowcolor=[]
        
        ##choosing min, median, max bursts
        
        segstart=dose_st+(dose_end-dose_st)//2
        
        segend=segstart+self.traceduration
        
        
        
        
        burstdict=channelbursts[(channelbursts['Timestamp [µs]']>=segstart) 
                                & (channelbursts['Timestamp [µs]']<=segend)]
        
        print(burstdict.shape, 'burstdictshape')
        
        p=figure(width=500, height=200, title=str(channelid),
                 tools="pan,wheel_zoom,box_zoom,reset,save,hover")

        
        
        
        
        for index, row in burstdict.iterrows():
            
            spk=row['Spikes']
            
            
            if iso==True:

                color=map_value_to_color(row['iso_forest_scores'],
                                         self.burstwithoutlierscolor_mapper)
            
            else:
                color='red'

            times=[int(time_reverse(x, dose_st)) for x in spk]
            
            boxleft=times[0]
            boxright=times[-1]
            
            boxesleft.append(boxleft)
            boxesright.append(boxright)
            shadowcolor.append(color)
            
            
            
        rawsxs=np.arange(int(time_reverse(segstart, dose_st)), int(time_reverse(segend, dose_st)), 1)
        
        rawsys=fv[int(time_reverse(segstart, dose_st)):int(time_reverse(segend, dose_st))]
        
        

        p.line(rawsxs, rawsys, color='grey')
        
        hline1 = Span(location=threshold, dimension='width',
                      line_color='red', line_width=2, line_dash='dotted')
                
        hline2 = Span(location=threshold*(-1), dimension='width', 
                      line_color='red', line_width=2, line_dash='dotted')
        
        p.add_layout(hline1)
        
        p.add_layout(hline2)
        
        p.yaxis.axis_label ='ADC'
        
        print(len(boxesleft), 'lenbox')
        
        for left, right, color in zip(boxesleft, boxesright,  shadowcolor):
            box=BoxAnnotation(left=left, right=right, fill_alpha=0.4, fill_color='blue')
            
            p.add_layout(box)
        
        ###each channel is a single plot for a well.####
       
        
        #color_bar = ColorBar(color_mapper=self.burstwithoutlierscolor_mapper,
                             #location=(0, 0))
        
        # Add the color bar to the plot
        #p.add_layout(color_bar, 'right')
        
        custom_locations =rawsxs[::40000].copy()  # Adjust the slicing based on your data
        # Ensure custom_locations are native Python integers
        custom_locations = custom_locations.astype(int)
        # Step 2: Set the x-axis ticker
        p.xaxis.ticker = FixedTicker(ticks=[int(tick) for tick in custom_locations]) 
        # Ensure they are Python integers

        # Step 3: Generate new label locations
        num_ticks = len(custom_locations)  # Number of custom tick locations
        new_locations = np.arange(2, 2 + num_ticks * 2, 2).astype(str)  # Create string labels for ticks
        
        # Step 4: Create a dictionary for custom labels
        custom_labels = dict(zip([int(tick) for tick in custom_locations], new_locations))  # Convert ticks to int for dict keys

        # Assign custom labels to x-axis
        p.xaxis.major_label_overrides = custom_labels


        
        
       
        self.burstplot_list.append(p)
        
        
        
    
                      
        
        
        
        
        
    
    
    
    
    
       
        
        
        
    
        
    
    
    
    
    
    
    
    


