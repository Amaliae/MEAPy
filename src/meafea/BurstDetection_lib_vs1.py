#!/usr/bin/env python
# coding: utf-8

# In[12]:


import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from sys import exit
from lxml import etree
from scipy import signal 
from sklearn import preprocessing
from scipy import optimize
import regex as re
import timeit
from scipy.signal import butter, lfilter
from scipy.signal import freqz
from lmfit.models import StepModel, LinearModel
import sys, importlib, os
import McsPy.McsData
import McsPy.McsCMOS
from McsPy import ureg, Q_
from numpy import trapz
from sklearn.neighbors import LocalOutlierFactor
import h5py
import pickle
from sklearn.decomposition import PCA
from scipy.stats import skew
from scipy.stats import kurtosis
get_ipython().run_line_magic('matplotlib', 'inline')
import itertools
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm # to build a LOWESS model
from scipy.interpolate import interp1d 
import os
import changefinder
import nbimporter

from pathlib import Path

from meafea.MEAFEA_functions_vs1 import *
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
from bokeh.models import LinearColorMapper, ColorBar
output_notebook()


# In[13]:


help(zoom_burst)


# In[2]:


def find_peak(count):
    
    """Function to extract peaks in ISI distribution of the spike trains"""
    
    t1=np.diff(np.sign(np.diff(count))) #counts 
    peak_ind=np.where(t1==(-2))[0]+1
    diff_peak_ind=np.diff(peak_ind)
    
    
    ####print(peak_ind, 'peaks')
    ####print( diff_peak_ind, 'diff peaks')
    
    
   
    
    if len(diff_peak_ind[diff_peak_ind<=2])>0:
        
        to_delete=[]
        
        ####print('while loop')

        compinds=np.where(np.diff(peak_ind)<=2)[0].tolist()
        
        for compind in compinds:
            
            p_ind1=peak_ind[compind]
            p_ind2=peak_ind[compind+1]
            
            if count[p_ind1]>count[p_ind2]:
                
                to_delete.append(compind+1)
            else: 
                to_delete.append(compind)
                
        peak_ind=np.delete(peak_ind, to_delete)
        
    return peak_ind


# In[3]:


def peak(count, bins, th):
    
    """Function to extract peaks in ISI distribution of the spike trains"""
    
    t1=np.diff(np.sign(np.diff(count))) #counts 
    peak_ind=np.where(t1==(-2))[0]+1
    diff_peak_ind=np.diff(peak_ind)
    
    
    ##print(peak_ind, 'peaks')
    ##print( diff_peak_ind, 'diff peaks')
    
    
   
    
    if len(diff_peak_ind[diff_peak_ind<=2])>0:
        
        to_delete=[]
        
        ##print('while loop')

        compinds=np.where(np.diff(peak_ind)<=2)[0].tolist()
        
        for compind in compinds:
            
            p_ind1=peak_ind[compind]
            p_ind2=peak_ind[compind+1]
            
            if count[p_ind1]>count[p_ind2]:
                
                to_delete.append(compind+1)
            else: 
                to_delete.append(compind)
                
        peak_ind=np.delete(peak_ind, to_delete)

            
            
            
            
            
    ###print(peak_ind, 'after filter')
            
    fp=0
    sp=1
    spi=0
    fpi=0
    ind_peak=None
    ISIth=None
    mp=max(bins)
    
    
    ##print(bins[peak_ind], 'isis')
    
    if len(peak_ind)>=2:
        
        ###print('two peaks ', bins[peak_ind])
        
        test_bin=bins[peak_ind]
        
        
        if len(test_bin[test_bin<np.log10(th)])>0:
            
            ####print('yes, smaller')
            
        

            for index, ind in enumerate(peak_ind):
                
                ###print(index, 'index here')

                if (bins[ind]<np.log10(th)) and (count[ind]>fp): 
                    
                    ###print(ind, 'index existed')#finds max peak left than 11!! not adapitive

                    fp=count[ind] #peak count 
                    fpi=ind
                    ind_peak=index
                    
                    
                    # peak index 
                    
                    ###print(peak_ind[(ind_peak+1):], 'peakindexesremained')

                    for i in peak_ind[(ind_peak+1):]:  #for peak following burst peak find minimim count? strange
                        gmin=min(count[fpi:i])
                        
                        ####print(gmin, 'gmin')##finding the void  between burst and following peak

                        geomin=np.sqrt((count[fpi])*(count[i])) 
                        
                        ###print(geomin, 'geomin')#geometric mean of peaks

                        #this works for some reason: 
                        #could have been, peak with smallest geometric mean 

                        #the principle is no spikes between bursts (as if all neurons where involved)
                        
                        ###print(1-(gmin/geomin), 'geo')


                        if (1-(1-(gmin/geomin)) <sp) and (1-(gmin/geomin)>0.35): #if where two peak having most rare ISIs in between whose are burst peak and pause
                            sp=1-(1-(gmin/geomin))
                            spi=i
                            ISIth=bins[np.argmin(count[fpi:spi])]
                        
                            

        else:
            for index, ind in enumerate(peak_ind):
                
                ###print(index, 'index here')
                 
                if (bins[ind]<mp): 
                    
                    ###print(ind, 'index existed')#finds max peak left than 11!! not adapitive

                    mp=bins[ind] #peak count 
                    fpi=ind
                    ind_peak=index
                    
                    
                    # peak index 
                    
                    ###print(peak_ind[(ind_peak+1):], 'peakindexesremained')

                    for i in peak_ind[(ind_peak+1):]:  #for peak following burst peak find minimim count? strange
                        gmin=min(count[fpi:i])
                        
                        ###print(gmin, 'gmin')##finding the void  between burst and following peak

                        geomin=np.sqrt((count[fpi])*(count[i])) 
                        
                       # ##print(geomin, 'geomin')#geometric mean of peaks

                        #this works for some reason: 
                        #could have been, peak with smallest geometric mean 

                        #the principle is no spikes between bursts (as if all neurons where involved)
                        
                        ###print(1-(gmin/geomin), 'geo')


                        if (1-(1-(gmin/geomin)) <sp) and (1-(gmin/geomin)>0.35): #if where two peak having most rare ISIs in between whose are burst peak and pause
                            sp=1-(1-(gmin/geomin))
                            spi=i
                            ISIth=bins[np.argmin(count[fpi:spi])]
                        
                           
                            
    return (ISIth)


# In[26]:


def LogISI(sp, f, th):
    
    """Function to  detect burst-like spike trains
    sp: spiketrain
    th: threshold"""
    
    Bursts=[]
    
    bind=[]
    merge=[]
    

    ISI=np.log10(np.diff(sp))
    
    #f=len(sp)/((sp[-1]-sp[0])/1000000)
    
    
    
    test_surprise=[]
    #log values of spikes
    
    count, bins=np.histogram(ISI, bins=np.arange(0, np.max(ISI)+0.2, 0.1))
    bins=bins[:-1]
    y_hat= sm.nonparametric.lowess(exog=bins, endog=count, frac=0.3)[:, 1] #
    
    
    #filter
    
    ISIth=peak(y_hat, bins, th)
    
    ##print(ISIth, 'ISIth')
    
    if ISIth!=None:
        
        ##print('checking this')
        
        
        
        extendFlag=1
        maxISI1 = np.log10(th)
        maxISI2=ISIth

    
    
    
        if ISIth<np.log10(th):
            
            ##print('Case 1, Peak based')

            maxISI1=ISIth

            extendFlag=0

            cross_ids=np.where(ISI<maxISI1)[0] 
            ##print(len(cross_ids))###(threshold/8))[0]

            diffcross_ids=np.diff(cross_ids).tolist()

            for i in range(len(diffcross_ids)):

                if diffcross_ids[i] < 2:

                    bind=bind+cross_ids.tolist()[i:i+2]

                else:

                    if len(bind)>1:
                        Bursts.append(sp[min(bind):(max(bind)+2)].tolist())

                    bind=[]
                    continue

        else:
            
            ##print('Case 2' , 'long peak')

            cross_ids=np.where(ISI<maxISI1)[0] 
            ###print(len(cross_ids))###(threshold/8))[0]

            diffcross_ids=np.diff(cross_ids).tolist()

            for i in range(len(diffcross_ids)):

                if diffcross_ids[i] < 2:

                    bind=bind+cross_ids.tolist()[i:i+2]

                else:

                    if len(bind)>5:

                        right=bind[-1]
                        left=bind[0]


                        for index in range(left, 0, -1):
                            if ISI[index]<maxISI1:
                                bind=[index]+bind
                            else:
                                break
                        for index in range(right, len(ISI), 1):
                            if ISI[index]<maxISI2:
                                bind=bind+[index]
                            else:
                                break

                        Bursts.append(sp[min(bind):(max(bind)+2)].tolist())

                        bind=[]

                        continue

                    else:
                        bind=[]
                        continue


                   
        ##print(len(Bursts))

        if extendFlag==1:

            if len(Bursts) > 0:       
                merge.append(Bursts[0])
                
               

                for B in Bursts[1:]:

                    if B[0]-merge[-1][-1] < maxISI2: #merge burst that are not separated with void///update///
                        merge[-1].extend(B)
                        
                    else:
                        merge.append(B)
                       
                        
                        
            Bursts=merge
 
    else:
        ##print('Case 0, no peak')
        maxISI1=np.log10(th)

        cross_ids=np.where(ISI<maxISI1)[0] 
        ###print(len(cross_ids))###(threshold/8))[0]

        diffcross_ids=np.diff(cross_ids).tolist()

        for i in range(len(diffcross_ids)):

            if diffcross_ids[i] < 2:

                bind=bind+cross_ids.tolist()[i:i+2]

            else:

                if len(bind)>5:
                    Bursts.append(sp[min(bind):(max(bind)+2)].tolist())

                bind=[]
                continue

        
                
    return (Bursts, np.ones(len(Bursts)))
            
            


# In[2]:


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
            
            


# In[18]:


def segment_b(stamp, Label, Type, Duration=None):
    
    
    
    ##print(Label, 'Label')#??????
    
    if Type=='Spike':
        
        start=int(stamp[str([str(Label)])]['start'])
        end=int(stamp[str([str(Label)])]['stop'])
        ##print('end', end)
        ##print('start', start )
        
    elif Type=='Burst':
        start=int(stamp[str([str(Label)])]['start'])
        end=int(stamp[str([str(Label)])]['stop'])
        ##print('end', end)
        ##print('start', start )
        
        
    if Duration:
        
        end=start+(Duration*1000000)
    
    return (start, end)


# In[28]:


def state_change(surprises, sp):
    
    #is surprise positive value 
    
    ##if change is more than 50% new state? 
    ##threshold for 3std
    
    surprises=surprises[1:]
    
    ##print(len(surprises), 'lengthofsurprise')
    sp=sp[2:]
    
    ##print(len(sp), 'lengthofspikes')
    
    threshold=np.std(surprises)*3
    
    highstatetrans=(surprises>threshold).astype('int')
    
    lowstatetrans=(surprises<threshold).astype('int')
    
    ##print(highstatetrans, np.diff(highstatetrans), 'initinnf')
    
    highstatest=np.where(np.diff(highstatetrans)==1)[0]
    highstateend=np.where(np.diff(highstatetrans)==-1)[0]
    
    ##print(highstatest, highstateend, '1ns,-1ns')
    
    hst=[]
    hend=[]
    
    if len(highstatest)>0 and len(highstateend)>0:
        
        if highstateend[0]<highstatest[0]:
            
            hst.append(0)
            hend.append(highstateend[0])
    
    if len(highstatest)>0:
    
        for stindex in highstatest:



            pair=-1
            
            


            if len(highstateend)>0:
                
                for endindex in highstateend:
                    if endindex>stindex: ###if ends when started

                        pair=endindex


                        hst.append(stindex)
                        hend.append(endindex)

                        break




            if pair==-1:  #if start does not end

                hst.append(stindex)
                hend.append(len(surprises)-1)
    else:
        if len(highstateend)>0:
            
            hst.append(0)
            hend.append(highstateend[0])
    
            

            
            
            
        
        
    ##print(hst, hend, 'hst, hend')
    
    
    
    
        
    hstlow=[0]+hend
    
    
    hendlow=hst+[len(surprises)-1]
                
                
    ##print(hstlow, hendlow, 'lwosfg')
    
    highspikesindex=[np.arange(hst[ind], hend[ind]+1).tolist() for ind in range(len(hst))]
    
    lowspikesindex=[np.arange(hstlow[ind], hendlow[ind]+1).tolist() for ind in range(len(hstlow))]
    
    highspikes=[sp[hst[ind]:hend[ind]+1] for ind in range(len(hst))]
    
    lowspikes=[sp[hstlow[ind]:hendlow[ind]+1] for ind in range(len(hstlow))]
    
    y_hat= sm.nonparametric.lowess(exog=np.arange(0, len(surprises)), endog=surprises, frac=0.1)[:, 1]
    
    #plt.plot(surprises, color='blue')
    
    
    
    #peaks_indexes=find_peak(y_hat)
    
    
    cf = changefinder.ChangeFinder()
    scores = [cf.update(p) for p in surprises]
    plt.plot(scores, color='red')
    
    plt.show()
    
    
    

    
    return highspikesindex, lowspikesindex, threshold, highspikes, lowspikes
    
    
    
    
    
    


# In[ ]:





# In[11]:


def find_peak_slope(b):
    
    st=b[0]
    
    isis=np.diff(b)
    
    minisiind=np.argmin(isis)+1
    
    sloperise=(1/(isis/1000000))/((b[minisiind]-b[0])/1000000)
    
    slopefall=(1/(isis/1000000))/((b[-1]-b[minisiind])/1000000)
    
    
    
    return (sloperise, slopefall)
    


# In[29]:


def dynsurprise(sp, ISI, f):
    
    
    SpikeGroups=[]
    BigSurprise=[]
    test_bind=[0, 1]
    
    for index in range(0, len(sp)-1):
        
        
        test_bind[1]=index
        
        if test_bind[1]-test_bind[0]>100:
            
            workingsp=sp[test_bind[0]:test_bind[1]+1]
            
            highr, lowr, th, highspikes, lowspikes=state_change(BigSurprise, workingsp)
            
            high=[i for el in highr for i in el]
            low=[i for el in lowr for i in el]
            
            SpikeGroups.extend(highspikes)
            SpikeGroups.extend(lowspikes)
            
            
            
            
            
            
            test_bind[0]=index
            
          
            
            
            
            plt.figure(figsize=(8, 10))
    
            plt.plot(BigSurprise[1:])
            y = np.interp(high, np.arange(0, len(BigSurprise[1:]), 1), BigSurprise[1:])
            w= np.interp(low, np.arange(0, len(BigSurprise[1:]), 1), BigSurprise[1:])
            
            plt.axhline(y=th)
            plt.plot(high, y, ls="", marker="*", ms=15,  color="crimson")
            plt.plot(low, w, ls="", marker="*", ms=15,  color="blue")
            
            
    
            plt.show()
        
            BigSurprise=[]
        
        
        
        
        surpr=Surprise(sp, test_bind, ISI, f)
        
        BigSurprise.append(surpr)
        
    
    plt.figure(figsize=(8, 10))
    
    plt.plot(BigSurprise)
    
    plt.show()
    
    return surpr
    
    
    
        
        
        
        
    
    
    
    


# In[30]:


###reattendng poisson 
##what does relatively long interval mean. 


# In[31]:


def PoissonSurprisetrue(sp, surpriseth, f):
    
    """Function to  detect burst-like spike trains"""
    
    Bursts=[]
    
    #container
   
   
    sp=np.sort(sp) #sorted spikes
    
    ISI=np.diff(sp) #ISIs 
    ISImean=np.mean(ISI)
    ISIhalfmean=ISImean/2
    ISIdoublemean=2*ISImean
    test_surprise=[]
    
    
    #f=len(sp)/((sp[-1]-sp[0])/1000000) #N/seconds###frequency 
    
    #esimench=dynsurprise(sp, ISI, f)
    
    
    
    ###print(f, 'Average FIring rate')
    start_indexes=np.where(ISI<ISIhalfmean)[0]
    
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
                    
                    ###stop adding if long interval and (or?) already 10 forward spikes failed, or

                    if ISI[j]<ISIdoublemean and fwd_failed<3:
                        test_bind[1]=j
                        

                        
                        p=Surprise(sp, test_bind, ISI, f)
                        
                       
                        
                        #if extension maximizes the surprise, extend

                        if p>p0:
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
                            #if removing from the start increases surprise, remove 
                            if p>p0:
                                bind[0]=i
                                p0=p
                            #else left is defined
                            else:
                                #if defined burst is higher than the given threshold
                                if p0>surpriseth:
                                    test_surprise.append(p0)
                                    
                                    Bursts.append(sp[bind[0]:bind[1]+1].tolist())

                                    bind=[-1, -1]
                                else:
                                    bind=[-1, -1]

                                break
                        break 
                          
    return (Bursts, test_surprise)
            
            


# In[32]:


def PoissonSurprise_test0(sp, surpriseth):
    
    """Function to  detect burst-like spike trains"""
    
    Bursts=[] #container
   
   
    sp=np.sort(sp) #sorted spikes
    plots=0
    ISI=np.diff(sp) #ISIs 
    ISImean=np.mean(ISI)
    ISIhalfmean=ISImean/2
    ISIdoublemean=2*ISImean
    
    test_surprise=[]
    
    test_Nofspikes=[]
    
    test_spikingrate=[]
    
    figevent, axesevent=plt.subplots(1, 1, figsize=(30, 6))
    axesevent.eventplot(sp, colors='grey', linewidths=0.8, alpha=0.6, linelengths=0.8)
    cm_subsection = np.arange(0, 55, 5)
    from  matplotlib.colors import Normalize
    cm_normed=Normalize(vmin=0, vmax=50, clip=False)
    from matplotlib import cm
    
   
    colors=[ cm.jet(cm_normed(x)) for x in cm_subsection ]
    
    ##print(cm_normed(x)) for x in cm_subsection]
    
    
    
    
    
    
    f=len(sp)/((sp[-1]-sp[0])/1000000) #N/seconds
    
    esimench=dynsurprise(sp, ISI, f)
    
    ###print(f, 'Average FIring rate')
    start_indexes=np.where(ISI<ISIhalfmean)[0]
    
    end_index=-1
    
    for i in range(len(start_indexes)-3):
        
        if  start_indexes[i]>end_index:
            bind=[-1, -1]

            if all(np.diff(start_indexes[i:i+4])==1):

                bind[0]=start_indexes[i]
                bind[1]=start_indexes[i+3] #or 4
                

                p0=Surprise(sp, bind, ISI, f)
                
                
                fwd_failed=0
                
                test_bind=bind

                for j in range(start_indexes[i+3]+1, len(ISI)):
                    
                   

                    if ISI[j]<ISIdoublemean and fwd_failed<10:
                        test_bind[1]=j
                        

                        
                        p=Surprise(sp, test_bind, ISI, f)
                        
                       
                        
                        #if extension maximizes the surprise, extend

                        if p>p0:
                            bind[1]=j
                            p0=p
                            
                        #if extension is not maximizing, bursts right index is defined
                        elif p==p0:  
                            fwd_failed=fwd_failed+1
                            end_index=bind[1]
                        else:
                            fwd_failed=fwd_failed+10
                            end_index=bind[1]
                            
                            
                          
                            #search for left index
                    else:
                        end_index=bind[1]
                        for i in range(bind[0]+1, bind[1]):
                            p=Surprise(sp, [i, bind[1]], ISI, f)
                            #if removing from the start increases surprise, remove 
                            if p>p0:
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
                                    
                                    
                                    if (p0>=0) and (p0 <=50):
                                        
                                        color_index=int(p0//5)
                                        
                                        color=colors[color_index]
                                        
                                        
                                       
                                        
                                        axesevent.eventplot(sp[bind[0]:bind[1]+1].tolist(), colors=color, linewidths=0.8, linelengths=0.8)
                                        plt.annotate(str(p0), xy=(sp[bind[0]:bind[1]+1].tolist()[0], 1.5), xytext=(sp[bind[0]:bind[1]+1].tolist()[0], 1.5))
                                    
                                    
                                   
                                        
                                        
                                    Bursts.append(sp[bind[0]:bind[1]+1].tolist())

                                    bind=[-1, -1]
                                else:
                                    bind=[-1, -1]

                                break
                        break 
    
    
    from  matplotlib.colors import ListedColormap
    
    
    
    cmap =ListedColormap(colors)
    
    

    
    plt.colorbar(cm.ScalarMappable(norm=cm_normed, cmap=cmap), ax=axesevent)#, ticks=np.arange(0, 55, 5).tolist())
    
    plt.show()
    
  
    
    
    plot_burst(test_surprise, test_Nofspikes, test_spikingrate)

           

    return (Bursts, test_surprise)
            
            


# In[33]:


def return_bokeh_colormap(name):
    cm = plt.cm.get_cmap(name)
    colormap = [rgb_to_hex(tuple((np.array(cm(x))*255).astype('int'))) for x in range(0,cm.N)]
    return colormap
def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb[0:3]


# In[34]:


def PoissonSurprise(sp, f, surpriseth):
    
    
    ##print(surpriseth, 'supriseth')
    
    """Function to  detect burst-like spike trains"""
    
    Bursts=[] #container
   
   
    sp=np.sort(sp) #sorted spikes
    plots=0
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
                    
                   

                    if ISI[j]<ISIdoublemean and fwd_failed<10:
                        test_bind[1]=j
                        

                        
                        p=Surprise(sp, test_bind, ISI, f)
                        
                       
                        
                        #if extension maximizes the surprise, extend

                        if p>=p0:
                            bind[1]=j
                            p0=p
                            
                        #if extension is not maximizing, bursts right index is defined
                        
                            
                       
                        else:
                            fwd_failed=fwd_failed+3
                            end_index=bind[1]
                            
                            
                          
                            #search for left index
                    else:
                        end_index=bind[1]
                        
                        
                        for i in range(bind[0]+1, bind[1]):
                            
                          
                            
                            p=Surprise(sp, [i, bind[1]], ISI, f)
                            #if removing from the start only increases surprise, remove 
                            if p>1.5*p0:
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
            
            


# In[34]:


def PoissonSurpriseGraph(sp, f, surpriseth, failfactor, startfactor):
    
    
    ##print(surpriseth, 'supriseth')
    
    """Function to  detect burst-like spike trains"""
    
    Bursts=[] #container
   
   
    sp=np.sort(sp) #sorted spikes
    plots=0
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
                    
                   

                    if ISI[j]<ISIdoublemean and fwd_failed<10:
                        test_bind[1]=j
                        

                        
                        p=Surprise(sp, test_bind, ISI, f)
                        
                       
                        
                        #if extension maximizes the surprise, extend

                        if p>=p0:
                            bind[1]=j
                            p0=p
                            
                        #if extension is not maximizing, bursts right index is defined
                        
                            
                       
                        else:
                            fwd_failed=fwd_failed+3
                            end_index=bind[1]
                            
                            
                          
                            #search for left index
                    else:
                        end_index=bind[1]
                        
                        
                        for i in range(bind[0]+1, bind[1]):
                            
                          
                            
                            p=Surprise(sp, [i, bind[1]], ISI, f)
                            #if removing from the start only increases surprise, remove 
                            if p>1.5*p0:
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
            
            


# In[1]:


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


# In[35]:


def plot_burst(surpise, nofspikes, rateofspikes):
    
    figs, axes=plt.subplots(3, 1, figsize=(8, 10))
    
    count, bins=np.histogram(surpise, np.arange(0, 50, 5))
    
    axes[0].hist(count, bins=bins[:-1], histtype='step')
    axes[0].set_title('Histrogram of Surprises')
    axes[0].set_xlabel('Surprise')
    axes[0].set_ylabel('Number of Trains')
    
    axes[1].scatter(surpise, nofspikes)
    axes[1].set_title('Scatter of surprisea and spikes ')
    axes[1].set_xlabel('Surprise')
    axes[1].set_ylabel('Number of Spikes')
    
     
    axes[2].scatter(surpise, rateofspikes)
    axes[2].set_title('Scatter of surprisea and spike rate ')
    axes[2].set_xlabel('Surprise')
    axes[2].set_ylabel('Spiking rate')
    figs.tight_layout()

    
   
    
    
    plt.show()
    
    
    
    
    
    
    
    
    
        
        
    
    


# In[36]:


#modify poisson ina way that you don't remove from the start


# In[37]:


from scipy.stats import poisson


# 
# repeated the process until either a relatively long
# interval was encountered or inclusion of one or
# more (up to 10) additional spikes failed to increase
# the Poisson surprise. (The length of such a relatively
# long interval had been chosen to be twice the
# buffer-wide average.)

# In[38]:


def Surprise(sp, test_bind, ISI, f, external=False):
    
    
    #sp is a list of event times in microseconds
    
   #indexes of events in list



    if external==False:

        T=(sp[(test_bind[1]+1)]-sp[test_bind[0]])/1000000

        ###print(T, 'Interval')
        ###print((test_bind[1]-test_bind[0])+1, 'lentoolarge?')

        rng=(test_bind[1]-test_bind[0])
        
    else:
        
        T=(test_bind[-1]-test_bind[0])/1000000
        rng=len(test_bind)
        
        
    try:
        p=poisson.sf(rng, f*T)#cumulative=sum([((f*T)**i)/np.math.factorial(i) for i in range(rng+1)])
        
        #p=1-((np.exp(-(f*T)))*cumulative)
    
    except:
        ##print('overflow')
        p=poisson.sf(rng, f*T)
        


    #if (p<0) or (p>100000000):
        
        ##print(p, 'probabilty')
    

    surprise=-np.log2(p)
    
    ###print(surprise, 'suprise')
    
    

    
    return (surprise)


# f=5
# n=6
# T=0.4
# 
# 

# cumulative=sum([((f*T)**i)/np.math.factorial(i) for i in range(n+1)])
# 
# 
# p=1-((np.exp(-(f*T)))*cumulative)
# 
# ##print(p, -np.log10(p))

# In[39]:


###testing, calibrating Poisson Surprise

###PLot Surprise to Number of bursts
### Plot Surpise to Number of spikes in burst
### 


# In[40]:


#parameterssets

plogisi={'set1':[150000], 'set2':[100000], 'set3': [200000]}
pmaxinterval={'set2': [150000, 250000, 2500000, 3], 
            }
psurprise={'set1':[7]}

params_set={}

params_set['logisi']=plogisi
params_set['poisson']=psurprise
params_set['maxinterval']=pmaxinterval


# In[41]:


#parameterssets

plogisi={'set1':[150000]}
pmaxinterval={'set1':[100000, 150000, 200000, 250000, 3]}
psurprise={'set1':[10]}

params_set={}

params_set['logisi']=plogisi
params_set['poisson']=psurprise
params_set['maxinterval']=pmaxinterval


# In[42]:


def burst_stat(spikes, methods, to, file, stamp, params_set):
    
    
    
    """Function to extract burst parameters and write in table, 
    Inputs
    spikes: spikes dataframe, from single experiment
    methods: is list methods, possible- ['logisi', 'poisson', 'maxinterval']
    
    Outputs
    
    Write files in a folder, bursts per experment"""
    
    
    #spikes['Dose Label']=spikes['Dose Label'].apply(lambda x : 0 if x=='Control' else int(x)) #int
    

    df=pd.DataFrame()
    
    
    chids=np.unique(spikes['Channel ID'].values)
    
    compound=np.unique(spikes['Compound ID'])[0] 
    #Compound ID in later versions. 
    experiment=np.unique(spikes['Experiment'])[0]
    
   
   
    
    for method in methods:
        
        
        df=pd.DataFrame()
        
        
        if method=='logisi':
            
            
            
            for setn in list(params_set[method].keys()):
            
                #args=[150000]

                args=params_set[method][setn]

                dataset=create_dataset(spikes, chids, compound, experiment, df, stamp, LogISI, method, False,  *args)
                out_path = Path(to) / f"{file[:-10]}LogISIBursts{setn}.pkl"
                dataset[1:].to_pickle(out_path)
                
        if method=='poisson':
            
            #args=[5]
            
            for setn in list(params_set[method].keys()):
            
                #args=[150000]

                args=params_set[method][setn]
                
                #print(method, args,  'method, args', df.head())
            
                dataset=create_dataset(spikes, chids, compound, experiment, df, stamp, PoissonSurprise, method, False,  *args)
                out_path = Path(to) / f"{file[:-10]}PoissonBursts{setn}.pkl"
                dataset[1:].to_pickle(out_path)
                
            
        if method=='maxinterval':
            
            #args=[200000, 300000, 200000, 3]
            for setn in list(params_set[method].keys()):
                
                #print(method, args,  'method, args', df.head())
            
                #args=[150000]

                args=params_set[method][setn]
            
            
            
                dataset=create_dataset(spikes, chids, compound, experiment, df, stamp, MaxInterval, method, False,  *args)
                out_path = Path(to) / f"{file[:-10]}MaxIntervalBursts{setn}.pkl"
                dataset[1:].to_pickle(out_path)
                

            
            
    return df          
                
                


# In[43]:


def create_dataset(spikes, channels, compound, experiment, df, stamp, function, method, plot,  *args):
    
    '''Helper function to create a dataset'''
    
  
    params_zoom=[100, 2000, 8, True, 40]
    
    ##print(method)
   
    wells=np.unique(spikes['Well ID'].values)
    
    ps=[figure(width=800, height=500, title=method+str(wells[r])+experiment) for r in range(len(wells))]
    
    
                           
    for indx, chid in enumerate(channels):
        
        
        colors=return_bokeh_colormap('viridis')
        colors_range=np.linspace(0, 200, 256)
        div=np.diff(colors_range)[0]
        color_mapper = LinearColorMapper(palette = colors, low = 0, high=200)
        
    
        
        dlabels=np.unique(spikes[spikes['Channel ID']==chid]['Dose Label'].values)
        #print(spikes[spikes['Channel ID']==chid].iloc[0, 3:7])
        
        #print(dlabels, 'dlabels')
        for dlabel in dlabels:
            
            st, end=segment_b(stamp, dlabel, 'Burst')
            
            length=end-st
            
            #print(length, 'length')



            #print(compound, 'Compound')

            #print(dlabel, 'Label')

            #print(chid, 'Channel')

            spiketrain=spikes[(spikes['Channel ID']==chid) & (spikes['Dose Label']==dlabel)]['Timestamp [µs]'].sort_values(ascending=True).values
            wiid=spikes[(spikes['Channel ID']==chid) & (spikes['Dose Label']==dlabel)]['Well ID'].values[0]

            chlab=spikes[(spikes['Channel ID']==chid) & (spikes['Dose Label']==dlabel)]['Channel Label'].values[0]

            pind=np.where(wells==wiid)[0][0]
            ##print(pind, wells, wiid)
            g=spikes[spikes['Channel ID']==chid]['Well ID'].values[0]
            ##print(g, 'Well ID')#Well ID

            p1=ps[pind]
            
            p1.rect(x=spiketrain,  y=np.ones(len(spiketrain))*chlab, width=2000, height=0.5, fill_color="grey", line_color="blue")

            f=len(spiketrain)/(length/1000000)
            
            bursts_outputs=function(spiketrain, f, *args)
            

            bursts=bursts_outputs[0]
            
            ##print(bursts, 'bursts')

            burstspikes=sum([len(bur) for bur in bursts])




            

            ISI=np.diff(spiketrain)


            for burst in bursts:
                
                

                if len(burst)>0:
                    
                    #print('we have burs')

                    s=Surprise(spiketrain, burst, ISI, f, external=True)
                    
                    
                    if plot==True:
                        
                        #print('plot is true', burst[0]+((burst[-1]-burst[0])/2), spiketrain[0])
                    
                    
                    
                   
                        
                        #color_index=int(s//div)

                        #color=colors[color_index]


                        p1.rect(x=burst[0]+((burst[-1]-burst[0])/2), y=(np.ones(1)*chlab)+0.3, 
                                width=burst[-1]-burst[0], height=0.05,
                                fill_color='red', line_color='red')
                        
                        ##print('Im showing')
                        
                        
                        
                        
                        ps[pind]=p1
                        
                        ##print('again')
                        
                    b=np.array(burst)
                    tempdf=pd.DataFrame(spikes[spikes['Channel ID']==chid].iloc[0, 3:7]).T
                    tempdf['Dose Label']=dlabel
                    tempdf['Compound ID']=compound
                    tempdf['Experiment']=experiment
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
                    
                    res=zoom_burst(b, params_zoom)
                    for r in ['Rise', 'Fall', 'Area']:
                        
                        tempdf[r]=res[r]
                    
                    
                    
                    #print(tempdf, df, 'both')
                          
                          
                          
                          


                    df=pd.concat([df, tempdf], axis=0)
                    
                    #print(df.shape, 'dfshape')
        

    from  matplotlib.colors import ListedColormap
    
    if plot==True:
        

        cb = ColorBar(color_mapper = color_mapper, location=(5,6))

        [p.add_layout(cb, 'right') for p in ps ]              
        [show(p) for p in ps] 
    
    return df 


# In[15]:


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
    


# In[16]:


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
    
    


# In[17]:


def modify_b(b):
    
    if isinstance(b, str)==True:
    
        b=b[1:-1].split(' ')

        b=np.array([ch for ch in b if ch[:1].isnumeric()==True]).astype('float')
    
    
    return b


# In[9]:


def experiment_burst_detection(folder_path, save_to, identifier1, identifier2, method_set, params_set):
    
    """Function to detect bursts of single experiments and save tabular data, here folder_path=save_to (CSV)"""
    
    filenames=load_raw(folder_path, identifier1, identifier2)
    stamp_dict = np.load(str(folder_path / "stamp.npy"), allow_pickle=True)[()]
    
    ##print(filenames)
   
    
    for file in filenames:
        
        positionh5=file.index('h5')+2

        stamp=stamp_dict[file[:positionh5]]
        
        print(stamp, 'stamp', file[:positionh5])
        
        ##print(stamp, 'stamp')
        
       
        
        spikes=pd.read_csv(str(folder_path/file))
        
        cols=spikes.columns
        
        coluname=[col for col in cols if 'Unnamed' in col]
        
        spikes.drop(coluname, axis=1, inplace=True)
        
        #print(spikes.shape, 'spikes_shape')
        df=burst_stat(spikes, method_set, save_to, file, stamp,  params_set)
        
        
        
            
        
    return  filenames


# In[47]:


def load_raw(path, identifier, identifier2):
    
    dirs = os.listdir(path)
    
    
    filename=[]
    
   
    
    
    
    for f in dirs:
        
        if re.search('h5', f) and re.search(identifier, f) and re.search(identifier2, f):
            
            filename.append(f)
            
    
            
            
    return filename


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




