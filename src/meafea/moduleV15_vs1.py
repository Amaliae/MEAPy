#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from scipy.signal import find_peaks, peak_prominences
get_ipython().run_line_magic('matplotlib', 'inline')
import itertools
from bokeh.io import output_notebook, show
from bokeh.plotting import figure
import ipyparams
from scipy.stats import median_abs_deviation

output_notebook()


# In[2]:


### 

###make exception for corrupted files. 


# In[3]:


def Fold_Folders(mainpath):
    
    
    hdf5folder_path=mainpath+'\\'+'hdf5'
    
    metadatapath=mainpath+'\\'+'hdf5'
    
    save_to=mainpath+'\\'+'spike_detection_csv'
    
    
    if 'spike_detection_csv' not in os.listdir(mainpath):
        
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
        
    if len(walkpaths)>0:

        metadatapath=[]
        folder_path=[]

        
        
        for wpth in walkpaths:
            
            dirswpth=os.listdir(wpth)
            
            mwsdirpth=[i for i  in dirswpth if 'mws' in i ]
            
            hdf5dirpth=[i for i  in dirswpth if i.endswith('h5')==True]
            
            print(hdf5dirpth, 'hdf5dirpth')
            
            if len(mwsdirpth)>0:
                
                metadatapth=wpth
                metadatapath.append(metadatapth)
                
                
            if len(hdf5dirpth)>0:
                
                fdpath=wpth
                folder_path.append(fdpath)
                

    

    
                
    else:
        print('no mwsin directory')
        
  
        
        
    return (folder_path, metadatapath, save_to)
 
    


# In[4]:


#Raw Load function
def load_raw_old(path, meta_path,  identifier, identifier2):

    '''Functon that is used idependently but applied here abundantly, so removing it
    For use call load_raw_ind'''

    dirs = os.listdir(path)
    
    print(dirs, 'dirs')


    dirs=[f for f in dirs if (re.search('h5', f) and re.search(identifier, f) and (identifier2 in f))==True]
    
    return dirs 


# In[5]:


#Raw Load function
def load_raw0(path, meta_path,  identifier, identifier2):

    dirs = os.listdir(path)
    
    print(dirs, 'dirs')


    filename=[]
    experiments={}


    for f in dirs:

        if (re.search('h5', f) and re.search(identifier, f) and (identifier2 in f)):

            filename.append(f)
            MetadataFilename=meta_path+'/'+f[:-7]+'.mws'

            Recordingfile=path+'/'+f
            recording=McsPy.McsData.RawData(Recordingfile)

            tree = etree.parse(MetadataFilename)
            
            ExperimentID=tree.xpath('//ExperimentID/text()')[0]
            
            CompoundID=recording.comment

            time=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]

            
            if ExperimentID not in experiments.keys():

                experiments[ExperimentID]={CompoundID:[]}
            else:
                experiments[ExperimentID].update({CompoundID:[]})

    
       
            Ahigh=int(tree.xpath('//FrequencymHz/text()')[0])/1000
            Alow=int(tree.xpath('//FrequencymHz/text()')[1])/1000
    
         
            experiments[ExperimentID][CompoundID].extend(time[:].tolist())








    return filename, experiments


# In[6]:


#Raw Load function
def load_raw(path, meta_path,  identifier, identifier2):

    dirs = os.listdir(path)
    
    print(dirs, 'dirs')


    filename=[]
    experiments={}


    for f in dirs:

        if (re.search('h5', f) and re.search(identifier, f) and (identifier2 in f)):

            filename.append(f)
            MetadataFilename=meta_path+'/'+f[:-7]+'.mws'

            Recordingfile=path+'/'+f
            recording=McsPy.McsData.RawData(Recordingfile)

            tree = etree.parse(MetadataFilename)
            
            ExperimentID=tree.xpath('//ExperimentID/text()')[0]
            
            CompoundID=recording.comment

            time=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]

            
            if ExperimentID not in experiments.keys():

                experiments[ExperimentID]={CompoundID:[]}
            else:
                experiments[ExperimentID].update({CompoundID:[]})

    
       
            
    
         
            experiments[ExperimentID][CompoundID].extend(time[:].tolist())








    return filename, experiments


# In[7]:


def led_info(Recordingfilename, Metadatafile):

    """Fuction to extract the LED stimulated well label, led color, led times and duration of the pulse"""

    R=McsPy.McsData.RawData(Recordingfilename)

    Info={}

    LED_T=[]
    LED_D=[]
    W=[]
    L=[]


    if len(list(R.recordings[0].event_streams.keys()))>3:



        for i in list(R.recordings[0].event_streams[3].event_entity.keys()):



            Wells=R.recordings[0].event_streams[3].event_entity[i].info.__dict__['info']['SourceChannelLabels']

            Label=R.recordings[0].event_streams[3].event_entity[i].info.__dict__['info']['Label']

            #print(R.recordings[0].event_streams[3].event_entity[i].info.__dict__['info']['EventID'], 'led')
            LED_T.append(R.recordings[0].event_streams[3].event_entity[i].data[0].tolist())
            LED_D.append(R.recordings[0].event_streams[3].event_entity[i].data[1].tolist())

            W.append(Wells)
            L.append(Label)



        Info['Wells']=W
        Info['Label']=L
        Info['Time']=LED_T
        Info['Duration']=LED_D


    return Info


# In[8]:


def Dilutions(tree, Compound):
    r=tree.getroot()
    dilutions={}
    for compound in r[1][0][1].getchildren():
        Series=[]
        for ser in compound[0][-1][0][1].getchildren():
            Series.append(ser.attrib['Dilution'].split(" ")[0])

        dilutions[str(compound[0][1].text)]=Series
        
        
    print('dilutions', dilutions)

    Labels=dilutions[Compound]

    return Labels


# In[9]:


def lab_book(tree):

    labbook={}

    root=tree.getroot()


    lablist=root[1][1].getchildren()

    for i in range(len(lablist)):


        labbook[root[1][1][i].tag]=root[1][1][i].text


    return labbook





def sd_extension(mydict, i,  APs, ExperimentID, Welllabel, WellID, Channellabel, Compound, dose_label, es):

    size=len(APs[1])

    print(len(APs[1]),'times')

    print(len(APs[2]),'amplitude')

    print(len(APs[3]), 'duration')




    mydict['Timestamp [µs]'].extend(APs[1])
    mydict['Peak Amplitudes'].extend(APs[2])
    mydict['Duration'].extend(APs[3])
    mydict['Experiment'].extend(np.repeat(ExperimentID, size))
    mydict['Well Label'].extend(np.repeat(Welllabel[i], size))
    mydict['Well ID'].extend(np.repeat(WellID[i], size))
    mydict['Channel ID'].extend(np.repeat(i, size))
    mydict['Channel Label'].extend(np.repeat(Channellabel[i], size))
    mydict['Dose Label'].extend(np.repeat(dose_label, size))
    mydict['ES_condition'].extend(np.repeat(es, size))
    mydict['Compound ID'].extend(np.repeat(Compound, size))
    mydict['Threshold'].extend(np.repeat(APs[6][0], size))
    
    mydict['BaselineM'].extend(np.repeat(APs[6][1], size))
    
   
 
    return mydict


# In[12]:


def AP_detection_lofnewest(to, lof, signal, partlength, dead_time, threshold, negative, positive,
                           fs, duration, learn, noisereduced, start=0, wavefeature='False'):
    
    
    
    
    ##PAYTSARS SUGGESTION; SPLIT INTO MANY PARTS TAKE MIN OF DITRIBUTION
    split=int(len(signal)//(2*partlength))
    
    if split>0:
        
        parted_signal=np.array_split(signal, split)
        
        stds=[np.std(subsignal[:]) for subsignal in parted_signal]
    else:
        
        stds=[np.std(signal[:])]
        
    
    
    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration
    
    


    pretrig=(fs/1000)*1
    postrig=(fs/1000)*2

    ##maximum spike duration
    
    maxdur=duration*1.5
    


    ##thresholds
    #stds=[np.std(subsignal[:]) for subsignal in parted_signal]
    
    stdused=np.std(signal[:])
    
    minstd=min(stds)
    meanstd=np.mean(stds)
    
    medianstd=np.median(stds)
    
    stdlist=[minstd, medianstd]
    
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
        
        
    if noisereduced==True:
        
        noisikner=PipNoise_detection(signal, 0.5, 5, True, True, 20000, 3)[1]
        
        print(noisikner, 'noisikner', len(noisikner))
        
        
        if len(noisikner)>0:
            
            #print(noisikner, 'sectime')
            
            #indexes, 0.5ms pre post noise removed as well
            noisik=np.array([np.arange(n-500, n+500, 1) for n in noisikner]).ravel()
            
            #print(noisik, 'noisik')
            mins=np.array([m for m in mins if m not in noisik ])
            
            print(mins, 'after noisik ')
            
    if len(mins)>1:


        # for testing will be >=1 but for the data >1



        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)+1] for i in mins])
        
        print(spks_wavef.shape, duration, 'shapedur', len(mins))


        amplitudes=np.median(spks_wavef[:, :], axis=1)-spks_wavef[:, int(spks_wavef.shape[1]/2)]
        
        if wavefeature=='True':

            features=np.apply_along_axis(wvf_features, 1, spks_wavef)

            pca_features=pca(spks_wavef, 4)
          
            X=np.concatenate((amplitudes.reshape(amplitudes.shape[0], 1), pca_features), axis=1)
            
            print(X.shape,  spks_wavef.shape, 'shapes')
            
            
        else:
            
            features=np.empty([0, 24])
            pca_features=np.empty([0, 4])
            X=np.empty([0, 5])




        
 
       

    else:
        X=np.empty([0, 5])
        mins=np.array([])
        spks_wavef=np.empty([0, int(duration+1)])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])
        stdlist=[[], []]
        
        
        
    

    fig, ax = plt.subplots()
    ax.plot(signal)
    
    if len(mins)>0:


        y = np.interp(mins, np.arange(0, len(signal), 1), signal)
        ax.plot(mins, y, ls="", marker="*", ms=15,  color="crimson", alpha=0.3)
        plt.axhline(y=(-threshold*medianstd))
        plt.axhline(y=(threshold*medianstd))
    #plt.savefig(to+'/'+str(len(mins))+'Spikes'+'.tif', format='tif')
    plt.show()



    
    
        



    return (spks_wavef, ((mins*50)+(start*50)), amplitudes, X[:, 0], lof, 
            np.concatenate((features, pca_features), axis=1), stdlist)


# In[13]:


def peak(cut, positive, negative):

    sign_diff=np.diff(np.sign(np.diff(cut)))
    
    
    
    
    
    


    ###print(cut,  np.diff(cut),  np.sign(np.diff(cut)), sign_diff, 'signdiffpeak')

    peak_idx=[]


    if negative==True:

        peak_indices=(np.where(sign_diff==2)[0]+1).tolist()
        if len(peak_indices)>0:
            #all negative peaks  peaks
            pk=np.argmin(cut[peak_indices])
            peak_idx=[peak_indices[pk]]

    if (positive==True) and (len(peak_idx)==0):

        peak_indices=(np.where(sign_diff==-2)[0]+1).tolist()
        if len(peak_indices)>0:
            #all negative peaks  peaks
            pk=np.argmax(cut[peak_indices])
            peak_idx=[peak_indices[pk]]

        #maximum negative peak actual index
    #np.where((sign_diff==2)|(sign_diff==-2))[0].tolist()
    #So, really it would be nice to have both positive and negative peaks' cutouts and find AP via negative peak.


    return (peak_idx, peak_indices)


# In[14]:


def peaknew(cut, positive, negative):

    sign_diff=np.diff(np.sign(np.diff(cut)))
    

    ###print(cut,  np.diff(cut),  np.sign(np.diff(cut)), sign_diff, 'signdiffpeak')

    positive_index=[]
    
    
    negative_index=[]


    if negative==True:

        negative_index=(np.where(sign_diff==2)[0]+1).tolist()
       

    if (positive==True):

        positive_index=(np.where(sign_diff==-2)[0]+1).tolist()
        #maximum negative peak actual index
    #np.where((sign_diff==2)|(sign_diff==-2))[0].tolist()
    #So, really it would be nice to have both positive and negative peaks' cutouts and find AP via negative peak.


    return (negative_index, positive_index)


# In[15]:


def peakslope(cut, positive, negative):

    sign_diff=np.diff(np.sign(np.diff(cut)))
    

    ###print(cut,  np.diff(cut),  np.sign(np.diff(cut)), sign_diff, 'signdiffpeak')

    positive_index=[]
    
    
    negative_index=[]


    if negative==True:

        negative_index=(np.where(sign_diff==2)[0]+1).tolist()
       

    if (positive==True):

        positive_index=(np.where(sign_diff==-2)[0]+1).tolist()
        #maximum negative peak actual index
    #np.where((sign_diff==2)|(sign_diff==-2))[0].tolist()
    #So, really it would be nice to have both positive and negative peaks' cutouts and find AP via negative peak.


    return (negative_index, positive_index)


# In[ ]:





# In[16]:


def wvf_features(wvf):

    features=np.zeros(24)

    #shape


    FD=np.diff(wvf)
    sign_diff=np.diff(np.sign(FD))
    zero_cross=(np.where((sign_diff==-2)|(sign_diff==2))[0])
    #zer_cross=(np.where((sign_diff==2))[0]+1)
    ##print('zero cross', zer_cross)

    AP_peak_cross=int(len(wvf)//2)-1 #AP peak index in FD (zero_crossing)


    try:
        peaks=np.where(np.diff(np.sign(np.diff(FD)))==-2)[0]+1
    except:
        peaks=np.array([])
    try:
        valleys=np.where(np.diff(np.sign(np.diff(FD)))==2)[0]+1
    except:
        valleys=np.array([])
    try:
        p1=np.sort(peaks[(peaks-AP_peak_cross)<0])[-1]
    except:
        p1=0
    try:
        p2=np.argmin(FD)
    except:
        p2=0
    ### valley of FD
    try:
        p3=AP_peak_cross
    except:
        p3=0
        ## second zero crossing (crossing )
    try:
        p4=np.sort(peaks[(peaks-AP_peak_cross)>0])[0]
    except:
        p4=0
    try:#np.argmax(FD) #peak of PD
        p5=np.sort(zero_cross[(zero_cross-AP_peak_cross)>0])[0]
    except:
        p5=0#zero_crossing after AP_zero_cross
    try:
        p6=np.sort(valleys[(valleys-AP_peak_cross)>0])[0]
    except:
        p6=0#valley after zero cross


    #fig, ax = plt.subplots()
    #ax.plot(FD)

    #y = np.interp([p1, p2, p3, p4, p5, p6], np.arange(0, len(FD), 1), FD)
    #ax.plot([p1, p2, p3, p4, p5, p6], y, ls="", marker="*", ms=15,  color="crimson")
    #plt.title('Feature guides')
    #plt.show()



    try:
        features[0]=p5-p1
    except:
        features[0]=0

    try:
        features[1]= FD[p4]-FD[p2]
    except:

        features[1]=0

    try:
        features[2]=FD[p6]-FD[p2]  #F3 correlation with reference waveform is missing
    except:
        features[2]=0

    try:
        features[3]=np.log(np.abs((FD[p4]-FD[p2])/(p4-p2)))

    except:
        features[3]=0

    try:
        features[4]=(FD[p6]-FD[p4])/(p6-p4)

    except:
        features[4]=0

    try:
        features[5]=np.log(np.abs((FD[p6]-FD[p2])/(p6-p2)))

    except:
        features[5]=0

    try:
        features[6]=np.sqrt(np.abs(np.mean(FD[:p2]))) #root mean square of pre-event amplitudes (added np/abs to avoid neg.)

    except:
        features[6]=0
    try:
        features[7]=((FD[p2]-FD[p1])/(p2-p1))/((FD[p3]-FD[p2])/(p3-p2))

    except:
        features[7]=0

    try:
        features[8]=((FD[p4]-FD[p3])/(p4-p3))/((FD[p5]-FD[p4])/(p5-p4))

    except:
        features[8]=0

    try:
        features[9]=FD[p2]/FD[p4]

    except:
        features[9]=0

    #phase based features

    features[10:16]=FD[[p1, p2, p3, p4, p5, p6]]
    try:
        features[16:19]=np.diff(FD)[[p1-1, p3-1, p5-1]] #second derivative of  zeros_crossings
    except:
        features[16:19]=0


    #distribution based features

    features[19]=np.percentile(FD, 75)-np.percentile(FD, 25)
    features[20]=np.percentile(np.diff(FD), 75)-np.percentile(np.diff(FD), 25)
    features[21]=skew(FD)
    features[22]=skew(np.diff(FD))
    features[23]=kurtosis(FD)




    return features


# In[17]:


def pca(wvfs, comps):

    if wvfs.shape[0]>10:


        pca=PCA(n_components=comps)
        wvfs_pca=pca.fit_transform(wvfs)
    else:
        wvfs_pca=np.zeros([wvfs.shape[0], 4])



    return wvfs_pca


# In[18]:


#Clusterization







# In[ ]:





# In[19]:


def AP_cutouts(v, dead_time, threshold, negative, positive, fs, duration, start=0):

    """The spike detection function, which takes channel
    raw data as an input and outputs spike parameters, here the spike cutout also beign extracted, amplitude is absolute"""


    dv=np.diff(v)#the derivative of the signal

    sign_diff=np.diff(np.sign(dv)) #the difference of the signs of derivative

     #parameter to  filter spikes by duration

    dead_time=(dead_time*fs)/1000  #parameter for refractory period consideration


    if negative==False:  #detect positive peaks


        peak_idx=(np.where(sign_diff==(-2))[0]+1).tolist() # idx where the sign difference is -2, idx+1 is peak
        threshold_crossings_idx=np.where(v>(threshold)*np.std(v[:]))[0] #threshold filtering, via std of 1s of data
        AP_idx=np.intersect1d(peak_idx, threshold_crossings_idx)# intersect peaks and thresholds

    if positive==False:   #detect negative peaks

        peak_idx=(np.where(sign_diff==(2))[0]+1).tolist()#
        threshold_crossings_idx=np.where(v<(-threshold)*np.std(v[:]))[0]
        AP_idx=np.intersect1d(peak_idx, threshold_crossings_idx)

    if (positive==True) and (negative==True): #detect both, positive and negative peaks

        peak_idx=(np.where((sign_diff==(-2)) | (sign_diff==(2)))[0]+1).tolist()

        threshold_crossings_idx=np.where((v>(threshold)*np.std(v[:])) | (v<(-threshold)*np.std(v[:])))[0]

        AP_idx=np.intersect1d(peak_idx, threshold_crossings_idx)

        dead_idx=(np.where(np.diff(AP_idx)<dead_time)[0]+1).tolist()
        AP_idx=np.delete(AP_idx, dead_idx)

    while np.any(np.diff(AP_idx)<dead_time):

        dead_idx=(np.where(np.diff(AP_idx)<dead_time)[0]+1).tolist()
        AP_idx=np.delete(AP_idx, dead_idx)



    spks_wavef=np.array([v[i-int(duration/2):i+int(duration/2)] for i in AP_idx])

    features=np.apply_along_axis(wvf_features, 1, spks_wavef)
    pca_features=pca(spks_wavef, 4)

    fig, ax = plt.subplots()
    ax.plot(signal)


    y = np.interp(AP_idx, np.arange(0, len(v), 1), v)
    ax.plot(AP_idx, y, ls="", marker="*", ms=15,  color="crimson")
    plt.axhline(y=(-threshold*np.std(v[:])))
    plt.savefig(to+'/'+str(len(AP_idx))+'Spikes'+'.tif', format='tif')
    plt.show()

    return (AP_idx[1:], ((np.array(AP_idx)*50)+(start*50))[1:].tolist(),  features[1:, :], pca_features[1:, :])


# In[20]:


def thresholdsplot(signal, threshold):
    
    
    plt.figure(figsize=(8, 6))
    
    plt.plot(stds)
    plt.hlines(y=stdused, xmin=0, xmax=180, colors='red', label='Before')
    plt.hlines(y=meanstd, xmin=0, xmax=180, colors='black', label='Mean', linestyles='dotted')
    plt.hlines(y=minstd, xmin=0, xmax=180, colors='green', label='Min')
    plt.hlines(y=medianstd, xmin=0, xmax=180, colors='yellow', label='Median')
    
    plt.legend()
    
    plt.title('stdsvalues')
   
    plt.show()


# In[ ]:





# In[21]:


def AP_detection_lofnew(to, lof, signal, dead_time, threshold, negative, positive, fs, duration, learn, start=0,
                        
                        wavefeature='False'):
    
    
    ##PAYTSARS SUGGESTION; SPLIT INTO MANY PARTS TAKE MIN OF DITRIBUTION
    split=int(len(signal)//(2*20000)) ###those 2 seconds windows
    parted_signal=np.array_split(signal, split)
    
    
    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration
    
    


    pretrig=(fs/1000)*1
    postrig=(fs/1000)*2

    ##maximum spike duration
    
    maxdur=duration*1.5
    


    ##thresholds
    stds= [np.std(subsignal[:]) for subsignal in parted_signal]
    
    stdused=np.std(signal[:])
    
    minstd=min(stds)
    meanstd=np.mean(stds)
    
    medianstd=np.median(stds)
    
    stdlist=[minstd, medianstd]
    
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
   

   
    if len(mins)>1:


        # for testing will be >=1 but for the data >1



        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)+1] for i in mins])
        
        print(spks_wavef.shape, duration, 'shapedur', len(mins))


        amplitudes=np.median(spks_wavef[:, :], axis=1)-spks_wavef[:, int(spks_wavef.shape[1]/2)]

        if wavefeature=='True':

            features=np.apply_along_axis(wvf_features, 1, spks_wavef)

            pca_features=pca(spks_wavef, 4)
          
            X=np.concatenate((amplitudes.reshape(amplitudes.shape[0], 1), pca_features), axis=1)
            
            
            if len(mins)>500 and learn==True:
                lof.fit(X)

        

        else:
            features=np.empty([0, 24])
            pca_features=np.empty([0, 4])
            X=np.empty([0, 5])
            

        

    
           
       
 
       

    else:
        X=np.empty([0, 5])
        mins=np.array([])
        spks_wavef=np.empty([0, int(duration+1)])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])
        stdlist=[[], []]

    fig, ax = plt.subplots()
    ax.plot(signal)


    y = np.interp(mins, np.arange(0, len(signal), 1), signal)
    ax.plot(mins, y, ls="", marker="*", ms=15,  color="crimson", alpha=0.3)
    plt.axhline(y=(-threshold*medianstd))
    plt.axhline(y=(threshold*medianstd))
    #plt.savefig(to+'/'+str(len(mins))+'Spikes'+'.tif', format='tif')
    plt.show()

    
    print(spks_wavef.shape, len(((mins*50)+(start*50))),  len(amplitudes), np.concatenate((features, pca_features), axis=1).shape)




    return (spks_wavef, ((mins*50)+(start*50)), amplitudes, amplitudes, lof, 
            np.concatenate((features, pca_features), axis=1), stdlist)


# In[117]:


def AP_detection_lofmin(to, lof, signal, dead_time, threshold, negative, positive, fs, duration, learn, start=0,
                        
                        wavefeature='False'):
    
    
    ##PAYTSARS SUGGESTION; SPLIT INTO MANY PARTS TAKE MIN OF DITRIBUTION
    split=int(len(signal)//(5000)) #250us windows
    parted_signal=np.array_split(signal, split)
    
    
    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration
    
    


    pretrig=(fs/1000)*1
    postrig=(fs/1000)*2

    ##maximum spike duration
    
    maxdur=duration*1.5
    


    ##thresholds
    stds= [np.std(subsignal[:]) for subsignal in parted_signal]
    
    print(len(stds), 'length of the parted signal check for the parititon')
    
    stdused=np.std(signal[:])
    
    minstd=np.percentile(stds, 25)
    meanstd=np.mean(stds)
    
    medianstd=np.median(stds)
    
    stdlist=[medianstd, minstd]
    
    if positive==False:
        
        negpeak_idxs, pospeak_idxs=peaknew(signal, positive, negative)
        
        
        
        threshold_crossings_idx=np.where(signal <(-threshold*minstd))[0]
        
        #print(negpeak_idxs, pospeak_idxs, threshold_crossings_idx)
        
        peak_idx=np.intersect1d(negpeak_idxs, threshold_crossings_idx)

    else:
        negpeak_idxs, pospeak_idxs=peaknew(signal, positive, negative)
        
        threshold_crossings_idx=np.where((signal>(threshold)*minstd) | (signal<(-threshold)*minstd))[0]
        
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
   

   
    if len(mins)>1:


        # for testing will be >=1 but for the data >1



        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)+1] for i in mins])
        
        print(spks_wavef.shape, duration, 'shapedur', len(mins))


        amplitudes=np.median(spks_wavef[:, :], axis=1)-spks_wavef[:, int(spks_wavef.shape[1]/2)]

        if wavefeature=='True':

            features=np.apply_along_axis(wvf_features, 1, spks_wavef)

            pca_features=pca(spks_wavef, 4)
          
            X=np.concatenate((amplitudes.reshape(amplitudes.shape[0], 1), pca_features), axis=1)
            
            
            if len(mins)>500 and learn==True:
                lof.fit(X)

        

        else:
            features=np.empty([0, 24])
            pca_features=np.empty([0, 4])
            X=np.empty([0, 5])
            

        

    
           
       
 
       

    else:
        X=np.empty([0, 5])
        mins=np.array([])
        spks_wavef=np.empty([0, int(duration+1)])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])
        stdlist=[[], []]

    fig, ax = plt.subplots()
    ax.plot(signal)


    y = np.interp(mins, np.arange(0, len(signal), 1), signal)
    ax.plot(mins, y, ls="", marker="*", ms=15,  color="crimson", alpha=0.3)
    plt.axhline(y=(-threshold*minstd))
    plt.axhline(y=(threshold*minstd))
    #plt.savefig(to+'/'+str(len(mins))+'Spikes'+'.tif', format='tif')
    plt.show()

    
    print(spks_wavef.shape, len(((mins*50)+(start*50))),  len(amplitudes), np.concatenate((features, pca_features), axis=1).shape)




    return (spks_wavef, ((mins*50)+(start*50)), amplitudes, amplitudes, lof, 
            np.concatenate((features, pca_features), axis=1), stdlist)


# In[118]:


def detect_base_threshold(signal):
    
    ##PAYTSARS SUGGESTION; SPLIT INTO MANY PARTS TAKE MIN OF DITRIBUTION
    split=int(len(signal)//(100000)) #5 second windows
    parted_signal=np.array_split(signal, split)
    
    
    stds= [median_abs_deviation(subsignal[:]) for subsignal in parted_signal]
    
    medianstd=np.median(stds)
    
    
    return  medianstd
        
        


# In[10]:


def AP_detection_lofmad(to, lof, signal, dead_time, threshold, negative, positive, fs, 
                        duration, expernaltresh,  learn, start=0, wavefeature='False'):
    
    
    
    
    
    ##PAYTSARS SUGGESTION; SPLIT INTO MANY PARTS TAKE MIN OF DITRIBUTION
    
    
    split=int(len(signal)//(100000))#5 second windows
    
    if split<1:
        split=1
        
    parted_signal=np.array_split(signal, split)

    
    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration
    
    


    pretrig=(fs/1000)*1
    postrig=(fs/1000)*2

    ##maximum spike duration
    
    maxdur=duration*1.5
    


    ##thresholds
    stds= [median_abs_deviation(subsignal[:]) for subsignal in parted_signal]
    
    stdused=np.std(signal[:])
    
    minstd=min(stds)
    meanstd=np.mean(stds)
    
    medianstd=np.median(stds)
    
    stdlist=[medianstd, stdused]
    
    
    if expernaltresh!=(-1):
        
        print(expernaltresh, 'working with expernal')
        
        basestd=expernaltresh
        
    else:
        basestd=medianstd
        
    print(basestd, stdused, 'thresholdbasleine')
        
        
        
    if positive==False:
        
        negpeak_idxs, pospeak_idxs=peaknew(signal, positive, negative)
        
        
        
        threshold_crossings_idx=np.where(signal <(-threshold*basestd))[0]
        
        #print(negpeak_idxs, pospeak_idxs, threshold_crossings_idx)
        
        peak_idx=np.intersect1d(negpeak_idxs, threshold_crossings_idx)

    else:
        
        negpeak_idxs, pospeak_idxs=peaknew(signal, positive, negative)
        
        threshold_crossings_idx=np.where((signal>(threshold)*basestd) | (signal<(-threshold)*basestd))[0]
        
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
   

   
    if len(mins)>1:


        # for testing will be >=1 but for the data >1



        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)+1] for i in mins])
        
        print(spks_wavef.shape, duration, 'shapedur', len(mins))


        amplitudes=np.median(spks_wavef[:, :], axis=1)-spks_wavef[:, int(spks_wavef.shape[1]/2)]

        if wavefeature=='True':

            features=np.apply_along_axis(wvf_features, 1, spks_wavef)

            pca_features=pca(spks_wavef, 4)
          
            X=np.concatenate((amplitudes.reshape(amplitudes.shape[0], 1), pca_features), axis=1)
            
            
            if len(mins)>500 and learn==True:
                lof.fit(X)

        

        else:
            features=np.empty([0, 24])
            pca_features=np.empty([0, 4])
            X=np.empty([0, 5])
            

        

    
           
       
 
       

    else:
        X=np.empty([0, 5])
        mins=np.array([])
        spks_wavef=np.empty([0, int(duration+1)])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])
        stdlist=[[], []]

    
    print(spks_wavef.shape, len(((mins*50)+(start*50))),  len(amplitudes), np.concatenate((features, pca_features), axis=1).shape)




    return (spks_wavef, ((mins*50)+(start*50)), amplitudes, amplitudes, lof, 
            np.concatenate((features, pca_features), axis=1), stdlist)


# In[6]:


def AP_detection_lof(to, lof, signal, dead_time, threshold, negative, positive, fs, duration, learn, start=0):

    if positive==False:
        r1=(signal <=(-threshold*np.std(signal[10000:]))).astype('int64')

    else:
         r1=((signal <=(-threshold*np.std(signal[:]))) | (signal >=(threshold*np.std(signal[:])))).astype('int64')


    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration

    buffer=20 #1ms to detect aps with positive threshold crossings but not negative
    
    
    r1=r1[np.where(r1==0)[0][0]:]





    diffr1=np.diff(r1)
    ###print(r1[:5000])
    rights=np.where(diffr1==(-1))[0]+1
    lefts=np.where(diffr1==1)[0]

    #ranges=np.where(diffr1==(-1))[0]-np.where(diffr1==1)[0]
    lefts=lefts[:len(rights)]
    rights=rights[:len(lefts)] #


    #in range min index + left original index

    #dead time_sufficient
    dead_idx=(np.where(np.diff(lefts)<dead_time)[0]+1).tolist()
    lefts=np.delete(lefts, dead_idx)
    rights=np.delete(rights, dead_idx)
    ###print(len(lefts))

    """removing spikes that violate dead time criteria also 50us peaks"""

    while np.any(np.diff(lefts)<dead_time) or np.any((rights-lefts)<1) or np.any((rights-lefts)>60):

        dead_idx=(np.where(np.diff(lefts)<dead_time)[0]+1).tolist()
        lefts=np.delete(lefts, dead_idx)
        rights=np.delete(rights, dead_idx)
    #no data
        to_delmin=np.where((rights-lefts)<1)[0].tolist()
        lefts=np.delete(lefts, to_delmin)
        rights=np.delete(rights, to_delmin)

        to_delmax=np.where((rights-lefts)>60)[0].tolist()
        lefts=np.delete(lefts, to_delmax)
        rights=np.delete(rights, to_delmax)

    #try:

    ##print(lefts, 'lefts')
    ##print(rights, 'rights')

    mins=np.array([peak(signal[(lefts[i]):(rights[i]+1)], positive, negative)[0][0]+(lefts[i]) for i in range(len(lefts))
                   if len(peak(signal[(lefts[i]):(rights[i]+1)], positive, negative)[0])>0]) #future peak search
    ind=[i for i in range(len(lefts)) if len(peak(signal[lefts[i]:(rights[i]+1)], positive, negative)[0])>0]
    lefts=lefts[ind]
    rights=rights[ind]#future peak search

    indexes=[i for i in range(len(mins)) if (((mins[i]-int(duration/2))>0) & ((mins[i]+int(duration/2))<len(signal)))]
    mins=mins[indexes]
    lefts=lefts[indexes]
    rights=rights[indexes]



    #except:
        #ValueError



        #mins=[]


    if len(mins)>1:

        divide=(mins-lefts)
        divide[divide==0]=1
        rise_vel=(np.abs(signal[mins]-signal[lefts]))/divide





        divide=(rights-mins)
        divide[divide==0]=1
        fall_vel=(np.abs(signal[rights]-signal[mins]))/divide

        # for testing will be >=1 but for the data >1



        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)] for i in mins])


        amplitudes=spks_wavef[:, 0]-spks_wavef[:, int(spks_wavef.shape[1]/2)]

        features=np.apply_along_axis(wvf_features, 1, spks_wavef)



        pca_features=pca(spks_wavef, 4)



        X=np.concatenate((rise_vel.reshape(rise_vel.shape[0], 1), fall_vel.reshape(rise_vel.shape[0], 1),
                          amplitudes.reshape(amplitudes.shape[0], 1)), axis=1)
        #strange_index=np.where(fall_vel>10)[0].tolist()
        #strange=mins[strange_index]
        ###print(mins)



        if len(mins)>500 and learn==True:
            lof.fit(X)
  #u  outliers
        ###print(mins)

        if len(mins)>0:

            plt.figure()
            plt.xlabel('rise_vel')
            plt.ylabel('Fall_vel')
            plt.scatter(rise_vel, fall_vel)
            plt.scatter(X[:, 0], X[:, 1], color='yellow', alpha=0.2)
            plt.show()


    else:
        X=np.empty([0, 2])
        mins=np.array([])
        spks_wavef=np.empty([0, 60])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])

    fig, ax = plt.subplots()
    ax.plot(signal)


    y = np.interp(mins, np.arange(0, len(signal), 1), signal)
    ax.plot(mins, y, ls="", marker="*", ms=15,  color="crimson")
    plt.axhline(y=(-threshold*np.std(signal[:])))
    #plt.savefig(to+'/'+str(len(mins))+'Spikes'+'.tif', format='tif')
    plt.show()


    #plt.figure()
    #plt.plot(spks_wavef[0, :])
    #plt.show()

    #plt.figure()
    #plt.hist(rise_vel, bins='auto')
    #plt.show()

    #plt.figure()
    #plt.xlabel('rise_vel')
    #plt.ylabel('Fall_vel')
    #plt.scatter(rise_vel, fall_vel)
    #plt.scatter(X[:, 0], X[:, 1], color='yellow', alpha=0.2)

    #plt.savefig(to+'/'+str(len(mins))+'Features'+'.tif', format='tif')

    #only returns rise_vel but fall _vell is also important
    #Future spike wvf feature extraction , classification
    ##print(features.shape,'rfshape' )

    ##print(pca_features.shape,'pcarfshape' )




    return (spks_wavef, ((mins*50)+(start*50)), amplitudes, X[:, 0], lof, np.concatenate((features, pca_features), axis=1), threshold*np.std(signal[:]))


# In[7]:


def AP_detection_lofNoise(signal, baseline, dead_time, threshold, negative, positive, fs, duration, start=0):

    if positive==False:
        r1=(signal <=(-threshold*np.std(signal[:]))).astype('int64')

    else:
         r1=((signal <=(-threshold*np.std(signal[:]))) | (signal >=(threshold*np.std(signal[:])))).astype('int64')


    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration

    pretrig=(fs/1000)*1
    postrig=(fs/1000)*2

    buffer=20 #1ms to detect aps with positive threshold crossings but not negative

    r1=r1[np.where(r1==0)[0][0]:]




    diffr1=np.diff(r1)

    ###pretrig is 2 data points, so there could be -1 or -2 at the the start of the trace if so replace with 0


    rights=(np.where(diffr1==(1))[0]+postrig).astype('int')
    lefts=(np.where(diffr1==1)[0]-pretrig).astype('int')

    lefts[(np.where(lefts)<0)[0]]=0



    ###

    true_rights=np.where(diffr1==(-1))[0]+1
    true_lefts=np.where(diffr1==1)[0]


    #ranges=np.where(diffr1==(-1))[0]-np.where(diffr1==1)[0]
    lefts=lefts[:len(rights)]
    rights=rights[:len(lefts)] #
    #print(lefts, rights, 'before')

    true_lefts=true_lefts[:len(true_rights)]
    true_rights=true_rights[:len(true_lefts)]

    #in range min index + left original index

    #dead time_sufficient
    dead_idx=(np.where(np.diff(lefts)<dead_time)[0]+1).tolist()

    dead_idx=(np.where(np.diff(true_lefts)<dead_time)[0]+1).tolist()

    #print(dead_idx)
    lefts=np.delete(lefts, dead_idx)
    rights=np.delete(rights, dead_idx)

    true_lefts=np.delete(true_lefts, dead_idx)
    true_rights=np.delete(true_rights, dead_idx)

    #print(lefts, rights, dead_time, rights-lefts, 'deadidx')
    ##print(len(lefts))

    """removing spikes that violate dead time criteria also 50us peaks"""

    while np.any(np.diff(lefts)<dead_time) or np.any((rights-lefts)<1) or np.any((rights-lefts)>maxdur):

        dead_idx=(np.where(np.diff(lefts)<dead_time)[0]+1).tolist()
        lefts=np.delete(lefts, dead_idx)
        rights=np.delete(rights, dead_idx)
    #no data
        to_delmin=np.where((rights-lefts)<1)[0].tolist()

        #print(to_delmin, 'delmin')
        lefts=np.delete(lefts, to_delmin)
        rights=np.delete(rights, to_delmin)

        to_delmax=np.where((rights-lefts)>maxdur)[0].tolist()

        #print(to_delmax, 'delmax')
        lefts=np.delete(lefts, to_delmax)
        rights=np.delete(rights, to_delmax)

    #try:

    #print(lefts, 'lefts', 'after')
    #print(rights, 'rights')

    mins=np.array([peak(signal[(lefts[i]):(rights[i])], positive, negative)[0][0]+(lefts[i]) for i in range(len(lefts))
                   if len(peak(signal[(lefts[i]):(rights[i])], positive, negative)[0])>0]) #future peak search
    ind=[i for i in range(len(lefts)) if len(peak(signal[lefts[i]:(rights[i]+1)], positive, negative)[0])>0]
    lefts=lefts[ind]
    rights=rights[ind]#future peak search


    #print(mins, 'mins')
    indexes=[i for i in range(len(mins)) if (((mins[i]-int(duration//2))>0) & ((mins[i]+int(duration//2))<len(signal)))]

    #print(mins, 'mins')


    mins=mins[indexes]
    lefts=lefts[indexes]
    rights=rights[indexes]



    #except:
        #ValueError



        #mins=[]


    if len(mins)>1:

        divide=(mins-true_lefts)
        divide[divide==0]=1
        rise_vel=(np.abs(signal[mins]-signal[true_lefts]))/divide





        divide=(true_rights-mins)
        divide[divide==0]=1
        fall_vel=(np.abs(signal[true_rights]-signal[mins]))/divide

        # for testing will be >=1 but for the data >1



        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)] for i in mins])


        amplitudes=spks_wavef[:, 0]-spks_wavef[:, int(spks_wavef.shape[1]/2)]

        features=np.apply_along_axis(wvf_features, 1, spks_wavef)



        pca_features=pca(spks_wavef, 4)



        X=np.concatenate((rise_vel.reshape(rise_vel.shape[0], 1), fall_vel.reshape(rise_vel.shape[0], 1),
                          amplitudes.reshape(amplitudes.shape[0], 1)), axis=1)


        lof = load_model('Novelty.pkl')
        yhat = lof.predict(X)



        ##find novel peaks indexes

        #mask = yhat != -1






        #X=X[mask, :]
        #spks_wavef=spks_wavef[mask, :]
        #amplitudes=amplitudes[mask]
        ##print(mins, '3rd')



        #mins=mins[mask]

        #print(mins)







  #u  outliers
        ##print(mins)

        if len(mins)>0:

            plt.figure()
            plt.xlabel('rise_vel')
            plt.ylabel('Fall_vel')
            plt.scatter(rise_vel, fall_vel)
            plt.scatter(X[:, 0], X[:, 1], color='yellow', alpha=0.2)
            plt.show()


    else:
        X=np.empty([0, 2])
        mins=np.array([])
        spks_wavef=np.empty([0, 60])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])

    fig, ax = plt.subplots()
    ax.plot(signal)


    y = np.interp(mins, np.arange(0, len(signal), 1), signal)
    ax.plot(mins, y, ls="", marker="*", ms=15,  color="crimson")
    plt.axhline(y=(-threshold*np.std(signal[:])))
    #plt.savefig(to+'/'+str(len(mins))+'Spikes'+'.tif', format='tif')

    #plt.xlim([19950, 20050])
    plt.show()




   


# In[8]:


def PipNoise_detection(signal, dead_time, threshold, negative, positive, fs, duration):

    if positive==False:
        r1=(signal <=(-threshold*np.std(signal[:]))).astype('int64')

    else:
         r1=((signal <=(-threshold*np.std(signal[:]))) | (signal >=(threshold*np.std(signal[:])))).astype('int64')


    dead_time=(fs/1000)*dead_time

    duration=(fs/1000)*duration

    buffer=20 #1ms to detect aps with positive threshold crossings but not negative




    diffr1=np.diff(r1)
    ##print(r1[:5000])
    rights=np.where(diffr1==(-1))[0]+1
    lefts=np.where(diffr1==1)[0]

    #ranges=np.where(diffr1==(-1))[0]-np.where(diffr1==1)[0]
    lefts=lefts[:len(rights)]
    rights=rights[:len(lefts)] #


    #in range min index + left original index

    #dead time_sufficient
    dead_idx=(np.where(np.diff(lefts)<dead_time)[0]+1).tolist()
    lefts=np.delete(lefts, dead_idx)
    rights=np.delete(rights, dead_idx)
    ##print(len(lefts))

    """removing spikes that violate dead time criteria also 50us peaks"""

    #try:

    #print(lefts, 'lefts')
    #print(rights, 'rights')

    mins=np.array([peak(signal[(lefts[i]):(rights[i]+1)], positive, negative)[0][0]+(lefts[i]) for i in range(len(lefts))
                   if len(peak(signal[(lefts[i]):(rights[i]+1)], positive, negative)[0])>0]) #future peak search
    ind=[i for i in range(len(lefts)) if len(peak(signal[lefts[i]:(rights[i]+1)], positive, negative)[0])>0]
    lefts=lefts[ind]
    rights=rights[ind]#future peak search

    indexes=[i for i in range(len(mins)) if (((mins[i]-int(duration/2))>0) & ((mins[i]+int(duration/2))<len(signal)))]
    mins=mins[indexes]
    lefts=lefts[indexes]
    rights=rights[indexes]



    if len(mins)>1:

        divide=(mins-lefts)
        divide[divide==0]=1
        rise_vel=(np.abs(signal[mins]-signal[lefts]))/divide





        divide=(rights-mins)
        divide[divide==0]=1
        fall_vel=(np.abs(signal[rights]-signal[mins]))/divide

        # for testing will be >=1 but for the data >1
        
        spks_wavef=np.array([signal[i-int(duration/2):i+int(duration/2)] for i in mins])
        amplitudes=spks_wavef[:, 0]-spks_wavef[:, int(spks_wavef.shape[1]/2)]



        X=np.concatenate((rise_vel.reshape(rise_vel.shape[0], 1), fall_vel.reshape(rise_vel.shape[0], 1),
                              amplitudes.reshape(amplitudes.shape[0], 1)), axis=1)
        
        
        
        lof = load_model('Novelty.pkl')
        yhat = lof.predict(X)



        ##find novel peaks indexes
        masknoise = yhat == -1

        minsnoise=mins[masknoise]
        
        
        if len(minsnoise)>0:



            
            spks_wav=np.array([signal[i-int(duration/2):i+int(duration/2)] for i in minsnoise])
            amplitudes=spks_wav[:, 0]-spks_wav[:, int(spks_wav.shape[1]/2)]
            
        else:
            X=np.empty([0, 2])
            mins=np.array([])
            spks_wav=np.empty([0, 60])
            amplitudes=np.array([])
            features=np.empty([0, 24])
            pca_features=np.empty([0, 4])

            minsnoise=np.array([])



       

    else:
        X=np.empty([0, 2])
        mins=np.array([])
        spks_wav=np.empty([0, 60])
        amplitudes=np.array([])
        features=np.empty([0, 24])
        pca_features=np.empty([0, 4])

        minsnoise=np.array([])

    fig, ax = plt.subplots()
    ax.plot(signal, color='blue')

    if len(mins)>0:
 
        if len(minsnoise)>0:


            y = np.interp(minsnoise, np.arange(0, len(signal), 1), signal)
            ax.plot(minsnoise, y, ls="", marker="*", ms=15,  color="crimson")
            plt.axhline(y=(-threshold*np.std(signal[:])))
            plt.title('PIpnoiseplot')
            #plt.savefig(to+'/'+str(len(mins))+'Spikes'+'.tif', format='tif')
            plt.show()





    return (spks_wav, minsnoise, amplitudes, X, threshold*np.std(signal[:]))


# In[123]:


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


# In[124]:


def highpass(v, cutoff, order, fs):

    nyq = 0.5 * 20000

    normal_cutoff = cutoff/nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)

    filtered = signal.filtfilt(b, a, v)

    return filtered


# In[125]:


def step_fit(t, response):
    step_mod=StepModel(form='linear', prefix='step_', nan_policy='omit') #step function
    #index=np.where(np.array(response)==np.array(response.min()))[0][0] #we dont need this
    ##print(index)
    scale=response.std()/response.mean()
    pars = step_mod.make_params(step_center=t.mean(),
                         step_amplitude=response.std(),
                         step_sigma=np.abs(scale)*1.5)  # parameters for the step function sigma=1.4 best parameter
    out = step_mod.fit(response, pars, x=t)

    return out.best_fit


def salpa(t, *p):
    return p[0]+(p[1])*(t)+(p[2])*((t)**2)+(p[3])*((t)**3)  #curve function


def salpa_fit(v, centroid):

    plt.figure()
    plt.plot(v)
    plt.title('Raw')
    plt.show()


    d=len(v)-(centroid)  #get upper range of

    t=np.arange(-int(centroid/2), int(centroid/2)).astype('int64')
    #print(t.shape)#centered on zero

    response=v[:len(t)]
    #print(response.shape)#

    y=step_fit(t, response).tolist()

    nc=np.arange(int(centroid), d).astype('int64') #get centroids


    for n in nc:

        t=np.arange(-int(centroid), int(centroid)).astype('int64') #arange input from -centroid to centroid, around zero

        response=v[int(n-centroid):int(n+centroid)]  #get signal window



        params, params_covariance = optimize.curve_fit(salpa, t, response, p0=np.ones(4))
        #optimize the function

        y_nc=salpa(t, *params)  #predict the central point value

        y.append(y_nc[int(centroid)]) #append only central point to list


    dv=np.array(v[: len(y)])-np.array(y)

    plt.figure()
    plt.plot(dv)
    plt.show()


    return dv


# In[126]:


def Dose(event_entities):

    starts=[]
    stops=[]

    #print(event_entities.keys())

    for e in list(event_entities.keys()):
        ##print('key', e)


        label_info=event_entities[e].info.__dict__['info']['Label']
        ##print('label', label_info)

        time=event_entities[e].data[0]

        ##print('time', time)




        if re.search('ControlStart', label_info) or re.search('TestStart', label_info):
            starts.extend(time)

        if re.search('ControlStop', label_info) or re.search('TestStop', label_info):
            stops.extend(time)

    return starts, stops


# In[ ]:





# In[127]:


def RArtifact(v, centroid, cutoff, order, fs, alg='StepSalpa'):

    #p = figure(title="Electrical Noise", x_axis_label='x', y_axis_label='y')



    if alg=='StepSalpa':


        dv=salpa_fit(v, centroid)

    else:


        dv=highpass(v, cutoff, order, fs)

    #p.line(np.arange(0, len(v), 1), v, legend_label="noise", line_width=2, color='red')
    #p.line(np.arange(0, len(dv), 1), dv, legend_label="Salpa", line_width=2, color='blue')
    #show(p)
    return dv


# In[ ]:





# In[9]:


def SpikeDetection_EStim(to, lof, recording, Compound,  MetadataMaskPath, tree, ExperimentID,  ES_Electrodes, 
                         experiment, Analyze,
                          
                    detmethod, E_Stim,
                         fs, cutoff, highpassth, order, dead_time, threshold,
                          usewelllist, negative, positive,
                          PipNoise, noallwells, mindur, maxdur, recordinthreshold, onlyLed, detsortfeatures, Dosethresholddoses=(-1, 1)):
  


    """Function to detect spike and write spike detections in a table for a single recording file"""

    D={'Timestamp [µs]':[], 'Peak Amplitudes':[], 'Duration':[], 'Channel ID':[], 'Well ID':[], 'Well Label':[],
       'Channel Label':[], 'Experiment':[], 'Dose Label':[], 'Compound ID':[], 'ES_condition':[], 
       'Threshold':[], 'BaselineM':[]}
    
    
    F=np.empty([0, 28])
    W=np.empty([1, 1])

    
    #recording=McsPy.McsData.RawData(Recordingfilename)
    #tree = etree.parse(MetadataFilename)
    WellID=np.load(MetadataMaskPath+'/'+'WellID.csv.npy')
    Welllabel=np.load(MetadataMaskPath+'/'+'WellLabel.csv.npy')
    Channellabel=np.load(MetadataMaskPath+'/'+'ChannelLabel.csv.npy')

    ChannelID=list(recording.recordings[0].analog_streams[0].channel_infos.keys())

    ChannelID=[int(i) for i in ChannelID]
    #Compounds=tree.xpath('//CompoundID/text()')
    #print(Compounds)
    #Compound=recording.comment

    Labels=Dilutions(tree, Compound)
    
    print(Labels, Compound,  'Labels', 'Compound')
    lab=lab_book(tree)

    lof_trained=0



    ES_Times=[]
    ES_Amplitudes=[]
    ES_Durations=[]

   
    ES_Electrodes=[int(i) for i in ES_Electrodes]
    ES_Wells=np.unique(np.array([WellID[int(i)] for i in ES_Electrodes]))
    ES_all=[np.where(WellID==i)[0].tolist()[:] for i in ES_Wells]
    ES_all=[item for sublist in ES_all for item in sublist]

    
    exp_times=sorted(experiment)
    if onlyLed==True:
        exp_times=np.unique(np.array(exp_times)).tolist()


    
    starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
    stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]
    
    
    scale=recording.recordings[0].analog_streams[1].channel_infos[0].adc_step.magnitude
    adcz=recording.recordings[0].analog_streams[1].channel_infos[0].get_field('ADZero')
    
    print(scale, 'scale', adcz, 'adcz')
    
    #signal=(recording.analog_streams[i].channel_data[:]-adcz)*scale
    
    #experiment start in us seconds
    esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
    print(esim, exp_times, 'esim, exp_times', ExperimentID, 'ExperimentID')
    #dose start in us
    
    

    if len(ES_Electrodes)>0:

        if len(list(recording.recordings[0].event_streams[0].event_entity.keys()))>0:

          
            
            ES_Times=find_ES(recording).astype('int64')

            ES_Amplitudes=tree.xpath('//Phase2Amplitude/text()')[0]

            ES_Amplitudes=np.array([int(i) for i in ES_Amplitudes.split(' ')])


            
            ES_imes=[(ES_Times[(ES_Times >=esim[i]) & (ES_Times<esim[i+1])]-esim[i])/50  if i!=len(esim)-1 else (ES_Times[ES_Times>esim[i]]-esim[i])/50  for i in range(len(esim))]
            ES_mes=[(ES_Times[(ES_Times >=esim[i]) & (ES_Times<esim[i+1])]) if i!=len(esim)-1 else 
                    (ES_Times[ES_Times>esim[i]])  for i in range(len(esim))]
            ES_Durations=tree.xpath('//Phase2Duration/text()')[0]
            ES_Durations=[int(i) for i in ES_Durations.split(' ')][0]


    dict_ES={}



    dict_ES['All_Electrodes']=ES_all  #all electrodes in stim wells
    dict_ES['Amplitudes']=ES_Amplitudes
    dict_ES['Wells']=ES_Wells
    dict_ES['Electrodes']=ES_Electrodes #only stim electrodes
    dict_ES['Times']=ES_Times

    stamp_dict={}

    sample_start=0
    
    print(starts, 'starts')

    for j in range(len(starts)):

        #dose_label=Labels[j]
        dose_start=int(starts[j])
        dose_stop=int(stops[j])

        dose_st=esim[j]
        #print('exp_times', exp_times)
        dose_label=[Labels[i] for i in range(len(exp_times)) if exp_times[i]==dose_st]
        
        print(dose_start, dose_stop, dose_label, 'dosesnew')
        dose_et=esim[j]+(50*(dose_stop-dose_start))
        stamp_dict[str(dose_label)]={'start':dose_st, 'stop':dose_et}
        #print(dose_start, dose_stop, dose_label)

        if len(ES_Times)>0:

            ES_T=(ES_imes[j]+sample_start).astype('int64')
            ES_es=(ES_mes[j]).astype('int64')



            sample_start=dose_stop

            if len(ES_T)>0:

                centroid=int((ES_T[-1]-ES_T[-2])/100)
                #print('centroid', centroid)

                window=int(ES_T[-1]-ES_T[-2])-centroid
                #print('window', window)

        artifact={}

        if Analyze=='True':

            for i in ChannelID:

                #print(i, 'channel id')
                
               ##Dosethresholddoses is tuple if first index value is not -1 detect adaptive coommon per channel
               ### if index 2 value is not -1 use single threshold for all
                
                if Dosethresholddoses[0]!=-1:
                    
                     
                    indthresh=int(recordinthreshold.recordings[0].analog_streams[0].channel_infos[i].row_index)

                    vthresh=recordinthreshold.recordings[0].analog_streams[0].channel_data[indthresh,
                                                                                    Dosethresholddoses[0]:Dosethresholddoses[1]]

                    fvthresh=highpass(vthresh, highpassth, 3, 20000)

                    Dosethreshold=detect_base_threshold(fvthresh)
                    
                elif  Dosethresholddoses[2]!=-1:
                    
                    Dosethreshold=Dosethresholddoses[2]
                    
                else:
                    
                    Dosethreshold=-1
                    
                    
                    
                index=int(recording.recordings[0].analog_streams[0].channel_infos[i].row_index)
                
                 
                #print(Welllabel[i], 'well')

                artifact[i]=[]




                if ((i in ES_all) and (len(ES_T)))>0:

                    if i in ES_Electrodes: ##future save this info
                        print('Stimulation')



                    v=recording.recordings[0].analog_streams[0].channel_data[index, dose_start:ES_T[0]]

                    fv=highpass(v, highpassth, 3, 20000)


                    APs=detmethod(to, lof, fv, dead_time, threshold, negative, positive,
                                              fs, mindur, Dosethreshold, False, start=dose_start, wavefeature=detsortfeatures)


                    D=sd_extension(D, i,  APs, ExperimentID, Welllabel, WellID, Channellabel, Compound, dose_label, 0)


                    v=recording.recordings[0].analog_streams[0].channel_data[index, (ES_T[-1]+window):dose_stop]



                    fv=highpass(v, highpassth, 3, 20000)




                    APs=detmethod(to, lof, fv, dead_time, threshold, negative, positive, fs, mindur, Dosethreshold,
                                  False, start=(ES_T[-1]+window))
                    
                        

                    D=sd_extension(D, i,  APs, ExperimentID, Welllabel, WellID, Channellabel, Compound, dose_label, 2)


                    for t in range(len(ES_T)):

                        if E_Stim==True:


                            time=ES_T[t]

                            t_l=ES_es[t]

                            if t!=len(ES_T)-1:

                                centroid=int((ES_T[t+1]-ES_T[t])/100)


                                window=int(ES_T[t+1]-ES_T[t])-centroid

                            else:
                                centroid=int((ES_T[t]-ES_T[t-1])/100)


                                window=int(ES_T[t]-ES_T[t-1])-centroid


                            centroid=50

                            v=recording.recordings[0].analog_streams[0].channel_data[index, (time):(time+window+centroid)]



                            dv=RArtifact(v, centroid, cutoff, order, fs, alg='StepSalpa')


                            #fv=highpass(dv, 5, 3, 20000)


                            maxart=max(abs(v))-min(abs(v))

                            artifact[i].append(maxart)



                            #APs=AP_detection_lof(to, lof, dv[100:220], dead_time,
                                                      #threshold, negative, positive, fs, 2, False, start=(time+100))



                            #D=sd_extension(D, i,  APs, ExperimentID, Welllabel, WellID, Channellabel, Compound, dose_label, t_l)

                            APs=detmethodEstim(to, [],  dv[100:], dead_time, threshold, negative,
                                                      positive, fs, 3, Dosethreshold, False, True, start=(time)+100)
                            
                            
                         

                            D=sd_extension(D, i,  APs, ExperimentID, 
                                           Welllabel, WellID, Channellabel, Compound, dose_label, t_l+11000)



                else:

                    #print('here')

                    print(dose_label, 'dose')

                    print(ExperimentID, 'experiment')
                    
                    print(dose_start, dose_stop, 'dose data points')
                    

                    v=recording.recordings[0].analog_streams[0].channel_data[index, dose_start:dose_stop]
                    #fv=butter_bandpass_filter(v, 10, 3500,  20000, order=3)
                    fv=highpass(v, highpassth, 3, 20000)
                    
                    mfv=np.mean(v)/np.mean(fv)
                    
                    print(mfv, 'ratio v to filtered')
                    
                    
                    APs=detmethod(to, lof, fv, dead_time, threshold,
                                         negative, positive, fs,  mindur, Dosethreshold,  True, start=dose_st/50, 
                                  wavefeature=detsortfeatures)

                    #print(dose_st/50, 'start_input')
                    #print(len(APs[1]))
                    
                    D=sd_extension(D, i,  APs, ExperimentID, Welllabel, WellID, Channellabel, Compound, dose_label, 0)
                    lof=APs[4]
                    features=APs[5]
                    #print(W.shape, 'W')
                    #print(APs[0].shape, 'New w shape')
                    if W.shape[1]>2:
                        W=np.concatenate((W, APs[0]), axis=0)
                    else:
                        W=APs[0]

                    F=np.concatenate((F, features), axis=0)






        lof_trained=lof

        dict_ES['Artifact Values']=artifact

        for key in list(D.keys()):
            
            print(len(D[key]), key)



    return (pd.DataFrame(D), dict_ES, lof_trained, lab, stamp_dict, F, W, {})


# In[ ]:





# In[130]:


#experiments=np.load(os.path.join(r'G:\AH\iN47\Pharm\spike_detection_csv', 'stamp.npy'), allow_pickle=True)[()]


# In[131]:


#doses=experiment_single_dose(folder_path, metadata_path, file, dose_label,  experiments)


# In[ ]:





# In[ ]:





# In[158]:


def experiment_single_dose(folder_path, metadata_path, file, dose_label,  experiments):
    
    MetadataFilename=metadata_path+'/'+file[:-7]+'.mws'
    
    Recordingfile=folder_path+'/'+file

    print(Recordingfile, 'RecordingFile', MetadataFilename, 'MetadataFilename')


    recording=McsPy.McsData.RawData(Recordingfile)


    tree = etree.parse(MetadataFilename)

    ExperimentID=tree.xpath('//ExperimentID/text()')[0]

    CompoundID=recording.comment
    
    experiment=experiments[ExperimentID][CompoundID]
    
    exp_times=sorted(experiment)
    
    Labels=Dilutions(tree, CompoundID)
    
   
    
    

    
    exp_times=np.unique(np.array(exp_times)).tolist()



    starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
    stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]

    #experiment start in us seconds
    esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
    print(exp_times, Labels,  esim)
    
    dose_st=[exp_times[i] for i in range(len(Labels)) if Labels[i]==dose_label][0]
    
    print(exp_times, Labels,  esim, dose_st, range(len(esim)))
    
    js=[index for index in range(len(esim)) if esim[int(index)]==dose_st]
    
    if len(js)==1:
        
        j=js[0]
    
        dose_start=starts[j]

        dose_end=stops[j]

        print(dose_label, dose_start, dose_end, 'doses threshold file and dose')
        
        
    else:
        
        'Print wrong doses'
        
        dose_start=0
        dose_end=0


    return (dose_start, dose_end)


# In[133]:


def find_ES(recording):

    ESstart_key=-1
    ESstop_key=-1

    for key in recording.recordings[0].event_streams[0].event_entity.keys():

        if recording.recordings[0].event_streams[0].event_entity[key].info.info['Label']=='Stimulator Activation Switch Start':

            ESstart_key=key

        if recording.recordings[0].event_streams[0].event_entity[key].info.info['Label']=='Stimulator Activation Switch Stop':

            ESstop_key=key






    if ESstart_key!=-1:

        print('Electrical stimulation existance')

        ES_times=recording.recordings[0].event_streams[0].event_entity[ESstart_key].data[0]

        ES_stoptimes=recording.recordings[0].event_streams[0].event_entity[ESstop_key].data[0]

        ES_duration=(ES_stoptimes-ES_times)[0]
        
        
    return ES_times


# In[134]:


def ES_info(recording, tree, experiments, WellID, Compound):
    
    ES_Times=[]
    Es_Time=[]
    ES_Amplitudes=[]
    ES_Durations=[]
    
    dict_ES={}
    
    ExperimentID=tree.xpath('//ExperimentID/text()')[0]
    
    ES_Electrodes=tree.xpath('//ChannelID/text()')
    ES_Electrodes=[int(i) for i in ES_Electrodes]
    ES_Wells=np.unique(np.array([WellID[int(i)] for i in ES_Electrodes]))
    ES_all=[np.where(WellID==i)[0].tolist()[:] for i in ES_Wells]
    ES_all=[item for sublist in ES_all for item in sublist]

    exp_times=sorted(experiments[ExperimentID][Compound])
    




    starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
    stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]
    esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
    #print(esim, exp_times)
    #dose start in us

    if len(ES_Electrodes)>0:

        if len(list(recording.recordings[0].event_streams[0].event_entity.keys()))>0:
            
            
            ES_Time=find_ES(recording).astype('int64')

            
            ES_Amplitudes=tree.xpath('//Phase2Amplitude/text()')[0]

            ES_Amplitudes=np.array([int(i) for i in ES_Amplitudes.split(' ')])


            ES_Times=ES_Time[0::1]
            ES_imes=[(ES_Times[(ES_Times >=esim[i]) & (ES_Times<esim[i+1])]-esim[i])/50  if i!=len(esim)-1 else (ES_Times[ES_Times>esim[i]]-esim[i])/50  for i in range(len(esim))]
            ES_mes=[(ES_Times[(ES_Times >=esim[i]) & (ES_Times<esim[i+1])]) if i!=len(esim)-1 else (ES_Times[ES_Times>esim[i]])  for i in range(len(esim))]
            ES_Durations=tree.xpath('//Phase2Duration/text()')[0]
            ES_Durations=[int(i) for i in ES_Durations.split(' ')][0]


    dict_ES['All_Electrodes']=ES_all  #all electrodes in stim wells
    dict_ES['Amplitudes']=ES_Amplitudes
    dict_ES['Wells']=ES_Wells
    dict_ES['Electrodes']=ES_Electrodes #only stim electrodes
    dict_ES['Times']=ES_Time
    
    
    return dict_ES


# In[135]:


def Raw(to, Recordingfilename, MetadataFilename, MetadataMaskPath, experiments, E_Stim, onlyLed):


    """Function to detect spike and write spike detections in a table for a single recording file"""

    
    
   

    recording=McsPy.McsData.RawData(Recordingfilename)

    tree = etree.parse(MetadataFilename)
    WellID=np.load(MetadataMaskPath+'/'+'WellID.csv.npy')
    Welllabel=np.load(MetadataMaskPath+'/'+'WellLabel.csv.npy')
    Channellabel=np.load(MetadataMaskPath+'/'+'ChannelLabel.csv.npy')

    ChannelID=list(recording.recordings[0].analog_streams[0].channel_infos.keys())

    ChannelID=[int(i) for i in ChannelID]
    Compounds=tree.xpath('//CompoundID/text()')
    #print(Compounds)
    Compound=recording.comment

    Labels=Dilutions(tree, Compound)
    lab=lab_book(tree)
    
    ExperimentID=tree.xpath('//ExperimentID/text()')[0]
    
    exp_times=sorted(experiments[ExperimentID][Compound])
    
    
    if onlyLed==True:
        exp_times=np.unique(np.array(exp_times)).tolist()



    starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
    stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]
    esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
    #print(esim, exp_times)
    #dose start in us
    
    dict_ES=ES_info(recording, tree, experiments)
    
    stamp_dict={}

    sample_start=0
    
    
    ##if drug experiment, list of dose starts and stops, so N(doses) arrays.

    for j in range(len(starts)):

        #dose_label=Labels[j]
        dose_start=int(starts[j])
        dose_stop=int(stops[j])

        dose_st=esim[j]
        #print('exp_times', exp_times)
        dose_label=[Labels[i] for i in range(len(exp_times)) if exp_times[i]==dose_st]
        dose_et=esim[j]+(50*(dose_stop-dose_start))
        stamp_dict[str(dose_label)]={'start':dose_st, 'stop':dose_et}
        #print(dose_start, dose_stop, dose_label)
       

        
        for i in ChannelID:
            
            index=int(recording.recordings[0].analog_streams[0].channel_infos[i].row_index)

              
            v=recording.recordings[0].analog_streams[0].channel_data[index, dose_start:dose_stop]
            
            ##indtead of yield could be spike detection function 
            
            
            yield  v 
                    

       


    


# In[136]:


def NSpikeDetection_EStim(to, lof, Recordingfilename, MetadataFilename, MetadataMaskPath, experiments, Analyze,
                          
                    detmethod, E_Stim,
                         fs, cutoff, highpassth, order, dead_time, threshold,
                          usewelllist, negative, positive,
                          PipNoise, noallwells, mindur, maxdur, onlyLed):


    """Function to detect spike and write spike detections in a table for a single recording file"""

    D={'Timestamp [µs]':[], 'Peak Amplitudes':[], 'Duration':[], 'Channel ID':[], 'Well ID':[], 'Well Label':[],
       'Channel Label':[], 'Experiment':[], 'Dose Label':[], 'Compound ID':[], 'ES_condition':[], 
       'Threshold':[], 'BaselineM':[]}
    
    
    F=np.empty([0, 28])
    W=np.empty([1, 1])

    recording=McsPy.McsData.RawData(Recordingfilename)

    tree = etree.parse(MetadataFilename)
    WellID=np.load(MetadataMaskPath+'/'+'WellID.csv.npy')
    Welllabel=np.load(MetadataMaskPath+'/'+'WellLabel.csv.npy')
    Channellabel=np.load(MetadataMaskPath+'/'+'ChannelLabel.csv.npy')

    ChannelID=list(recording.recordings[0].analog_streams[0].channel_infos.keys())

    ChannelID=[int(i) for i in ChannelID]
    
    if noallwells==True:
        
        print('is it true?')
        
        analysis_chlist=np.where(np.isin(WellID.astype('int'), usewelllist)==True)[0]
        
    else:
        
        analysis_chlist=ChannelID
        
        
    Compounds=tree.xpath('//CompoundID/text()')
    #print(Compounds)
    Compound=recording.comment

    Labels=Dilutions(tree, Compound)
    
    
    lab=lab_book(tree)

    lof_trained=0



    ES_Times=[]
    ES_Amplitudes=[]
    ES_Durations=[]

    ExperimentID=tree.xpath('//ExperimentID/text()')[0]
    ES_Electrodes=tree.xpath('//ChannelID/text()')
    ES_Electrodes=[int(i) for i in ES_Electrodes]
    ES_Wells=np.unique(np.array([WellID[int(i)] for i in ES_Electrodes]))
    ES_all=[np.where(WellID==i)[0].tolist()[:] for i in ES_Wells]
    ES_all=[item for sublist in ES_all for item in sublist]

    exp_times=sorted(experiments[ExperimentID][Compound])
    if onlyLed==True:
        exp_times=np.unique(np.array(exp_times)).tolist()



    starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1]
    stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2]
    esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
    #print(esim, exp_times)
    #dose start in us

    if len(ES_Electrodes)>0:

        if len(list(recording.recordings[0].event_streams[0].event_entity.keys()))>0:

            ES=recording.recordings[0].event_streams[0].event_entity

            key=list(ES.keys())[0]
            ES_Times=((ES[key].__dict__['data'][0])).astype('int64')

            ES_Amplitudes=tree.xpath('//Phase2Amplitude/text()')[0]

            ES_Amplitudes=np.array([int(i) for i in ES_Amplitudes.split(' ')])


            ES_Times=ES_Times[0::1]
            ES_imes=[(ES_Times[(ES_Times >=esim[i]) & (ES_Times<esim[i+1])]-esim[i])/50  if i!=len(esim)-1 else (ES_Times[ES_Times>esim[i]]-esim[i])/50  for i in range(len(esim))]
            ES_mes=[(ES_Times[(ES_Times >=esim[i]) & (ES_Times<esim[i+1])]) if i!=len(esim)-1 else (ES_Times[ES_Times>esim[i]])  for i in range(len(esim))]
            ES_Durations=tree.xpath('//Phase2Duration/text()')[0]
            ES_Durations=[int(i) for i in ES_Durations.split(' ')][0]


    dict_ES={}



    dict_ES['All_Electrodes']=ES_all  #all electrodes in stim wells
    dict_ES['Amplitudes']=ES_Amplitudes
    dict_ES['Wells']=ES_Wells
    dict_ES['Electrodes']=ES_Electrodes #only stim electrodes
    dict_ES['Times']=ES_Times


   

    stamp_dict={}

    sample_start=0
    
    artifacts={}

    for j in range(len(starts)):

        #dose_label=Labels[j]
        dose_start=int(starts[j])
        dose_stop=int(stops[j])

        dose_st=esim[j]
        #print('exp_times', exp_times)
        dose_label=[Labels[i] for i in range(len(exp_times)) if exp_times[i]==dose_st]
        dose_et=esim[j]+(50*(dose_stop-dose_start))
        stamp_dict[str(dose_label)]={'start':dose_st, 'stop':dose_et}
        #print(dose_start, dose_stop, dose_label)

        if len(ES_Times)>0:

            ES_T=(ES_imes[j]+sample_start).astype('int64')
            ES_es=(ES_mes[j]).astype('int64')



            sample_start=dose_stop

            if len(ES_T)>0:

                centroid=int((ES_T[-1]-ES_T[-2])/100)
                #print('centroid', centroid)

                window=int(ES_T[-1]-ES_T[-2])-centroid
                #print('window', window)

        artifact={}

        if Analyze=='True':
            
            
            
            #print(analysis_chlist, 'analysislist')

            for i in ChannelID:
                

                
                
                if i in analysis_chlist:


                    index=int(recording.recordings[0].analog_streams[0].channel_infos[i].row_index)

                    print(Welllabel[i], 'well')

                    artifact[i]=[]
                    
                    if PipNoise==True:
                        
                        print('here', index)

                        
                        v=recording.recordings[0].analog_streams[0].channel_data[index, dose_start:dose_stop]
                        
                        plt.figure()
                        plt.plot(v)
                        plt.title('raw')
                        plt.show()
                        #fv=butter_bandpass_filter(v, 10, 3500,  20000, order=3)
                        fv=highpass(v, highpassth, 3, 20000)

                        output=PipNoise_detection(v, 0.5, 5, True, True, 20000, 3)
                        noiseindexes=sorted(output[1])
                        noiseindexes=np.array(noiseindexes)
                        
                        
                       
                        
                        if len(noiseindexes)>0:
                            
                            wind=30000
                            
                            
                            
                            noiseindexesu=[noiseindexes[0]]+noiseindexes[np.where(np.diff(noiseindexes)>wind)[0]+1].tolist()
                            
                            artifact[i]=noiseindexesu
                            
                            print(noiseindexesu, 'noiseindexesu')
                            
                            
                            
                            for ni, n in enumerate(noiseindexesu):
                                
                                
                                
                                    
                                    
                                if ni!=len(noiseindexesu)-1:
                                    
                                    endn=noiseindexesu[ni+1]-wind
                                    
                                else:
                                    endn=len(v)
                                
                                if ni!=0:
                                    
                                    startn=noiseindexesu[ni-1]+wind
                                    
                                else:
                                    startn=0
                                    
                                print('heretoo')

                                

                                noise=v[n-wind:n+wind]
                                
                                signal1=v[startn:n-wind]
                                signal2=v[n+wind:endn]
                                
                                
                                plt.figure()
                                
                                plt.plot(noise)
                                
                                plt.title('noise')
                                plt.show()
                                
                                
                                startsig1=startn
                                startnoise=n-wind
                                startsig2=n+wind
                                
                                if len(signal1)<500:
                                    noise=np.concatenate([signal1, noise])
                                    signal1=[]

                                    startnoise=startsig1
                                    
                                if len(signal2)<500:
                                    

                                    noise=np.concatenate([noise, signal2])

                                    signal2=[]

                                    #startnoise=startsig2
                                    
                                    
                                
                                for sig, startnf  in zip([signal1, noise, signal2], [startsig1, startnoise, startsig2]):
                                    
                                    
                                    if len(sig)>0:
                                        
                                        sigclean=RArtifact(sig, 5, 100, 9, 20000, alg='esim')
                                        
                                        APs=AP_detection_lofnewest(to, lof, sigclean, 1000,  0.5, 4.8, True,True, 20000, 
                                                            3, False, True, start=dose_st/50+startnf)
                                        
                                        D=sd_extension(D, i,  APs, ExperimentID,
                                                   Welllabel, WellID, Channellabel, Compound, dose_label, 0)
                                        lof=APs[4]
                                        features=APs[5]

                                        if W.shape[1]>2:
                                            W=np.concatenate((W, APs[0]), axis=0)
                                        else:
                                            W=APs[0]

                                        F=np.concatenate((F, features), axis=0)
                                        
                                        

                   
                        
                        
    






    artifacts[str(dose_label)]=artifact
    
    
    lof_trained=lof

    dict_ES['Artifact Values']=artifact

    for key in list(D.keys()):

        print(len(D[key]), key)



    return (pd.DataFrame(D), dict_ES, lof_trained, lab, stamp_dict, F, W, artifacts)


# In[137]:


def save_model(clf, filename):
    with open(filename, 'wb') as f:
        pickle.dump(clf, f)


# In[138]:


def load_model(filename):
    with open(filename, 'rb') as f:
        clf = pickle.load(f)
    return clf


# In[139]:


def Dose_Label(experiments, folder_path, filenames):

    MetadataFilename=folder_path+'/'+file[:-7]+'.mws'

    Recordingfile=folder_path+'/'+file


# In[140]:


def notebook_name():



    currentNotebook = ipyparams.notebook_name

    return currentNotebook



# In[168]:


def experiment_spike_detection(folder, save_to, identifier, identifier2, identifier3,  MetadataMaskPath, threshold_dict,  Analyze,
                               detmethod=AP_detection_lofnew,  onlyLed=True,
                               E_Stim=False, fs=20000, cutoff=100, highpassth=5, order=9, dead_time=0.5, 
                               threshold=4.5, usewelllist=[], 
                               negative=True,
                               positive=True, PipNoise=False, noallwells=False, mindur=3, maxdur=500, 
                               detsortfeatures='False'):
    
    
    

    """Function to run spike detection on experimental folder and save detectios"""
    
    
    

    folder_paths, meta_path, save=Fold_Folders(folder)
    
    
    print(folder_paths, meta_path, save, 'new_paths')
    
    
  
        
   


    #filenames, experiments=load_raw(folder_path, meta_path,  identifier1, identifier1)
    ES_entry={}
    LED_entry={}
    LabBook={}
    stamp={}
    
    PipArtifacts={}
    
    experiments={}



    lof = LocalOutlierFactor(novelty=True, contamination=0.1)
    
    
    for folder_path in folder_paths:
        

        filenames=os.listdir(folder_path)
        
        print(folder_path, 'folder_path')
        
        if identifier3 not in folder_path:



            for file in filenames:

                if (re.search('mwd', file) and re.search(identifier, file)):



                    metafname=file[:-7]+'.mws'
                    
                    

                    meta_pathf=[mp for mp in meta_path if metafname in os.listdir(mp)]


                    print(meta_pathf, 'metasubfolder')

                    if len(meta_pathf)>0:


                        MetadataFilename=meta_pathf[0]+'/'+file[:-7]+'.mws'



                        Recordingfile=folder_path+'/'+file

                        print(Recordingfile, 'RecordingFile', MetadataFilename, 'MetadataFilename')


                        recording=McsPy.McsData.RawData(Recordingfile)


                        tree = etree.parse(MetadataFilename)

                        ExperimentID=tree.xpath('//ExperimentID/text()')[0]

                        CompoundID=recording.comment

                        print('CompoundID')
                        
                        
                         
                        try: 
                            time=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
                            
                        except AssertionError:
                            
                            
                            
                            print(f'Corrupted file, proceeding to next{Recordingfile}')
                            continue 
                        


                       
                        print(experiments, 'experiments in a loop', experiments.keys())



                        if ExperimentID not in experiments.keys():

                            experiments[ExperimentID]={CompoundID:[]}
                            experiments[ExperimentID][CompoundID].extend(time[:].tolist())


                        else:

                            if CompoundID not in experiments[ExperimentID].keys():

                                experiments[ExperimentID].update({CompoundID:[]})


                                experiments[ExperimentID][CompoundID].extend(time[:].tolist())

                            else:

                                experiments[ExperimentID][CompoundID].extend(time[:].tolist())



                    else:

                        print('no mws for i given file in main path')






    print(experiments, 'experiments')
        #print(filenames)

        #filenames, experiments1=load_raw(folder_path, meta_path, identifier1, identifier2)
        
    for folder_path in folder_paths:
        
        if identifier3 not in folder_path:
        
            filenames=os.listdir(folder_path)
            
            

            print(folder_path, 'folder_path')
            
            print(experiments, 'Experimentsdict')
            
            if 'File' in list(threshold_dict.keys()):

                thresholdfilename=threshold_dict['File']
                
                thresholdmetadatafilename=thresholdfilename[:-7]+'.mws'
                
                metadata_paththreshold=[mp for mp in meta_path if thresholdmetadatafilename in os.listdir(mp)][0]
                
                dose_label=threshold_dict['Dose Label']
                
                thresholddoses=experiment_single_dose(folder_path, metadata_paththreshold, thresholdfilename, dose_label,  
                                                      experiments)
                
                
            elif 'User Threshold' in list(threshold_dict.keys()):
                
                
                thresholddoses=(-1, -1,  threshold_dict['User Threshold'])
                
            else:
                
                thresholddoses=(-1, -1, -1)
                
                
                
                
                
                
            for file in filenames:



                if (re.search('mwd', file) and re.search(identifier, file)):

                    metafname=file[:-7]+'.mws'

                    print(metafname, 'metafname')


                    meta_pathf=[mp for mp in meta_path if metafname in os.listdir(mp)]

                    print(file, meta_pathf, 'files')

                    if len(meta_pathf)>0:


                        MetadataFilename=meta_pathf[0]+'/'+file[:-7]+'.mws'


                        Recordingfile=folder_path+'/'+file
                        
                        recording=McsPy.McsData.RawData(Recordingfile)
                        
                        try: 
                            time=recording.recordings[0].analog_streams[1].timestamp_index[:, 0]
                            
                        except AssertionError:
                            
                            
                            print(f'Corrupted file, proceeding to next{Recordingfile}')
                            continue 
                        
                       
                        tree = etree.parse(MetadataFilename)


                        ExperimentID=tree.xpath('//ExperimentID/text()')[0]

                        CompoundID=recording.comment
                        print(CompoundID, MetadataFilename, 'CompoundID, meta')


                        Recordingfile=folder_path+'/'+file

                        ES_Electrodes=tree.xpath('//ChannelID/text()')

                        Ahigh=int(tree.xpath('//FrequencymHz/text()')[0])/1000
                        Alow=int(tree.xpath('//FrequencymHz/text()')[1])/1000


                        LED_entry[file]=led_info(Recordingfile, MetadataFilename)
                        
                        
                         
                        experiment=experiments[ExperimentID][CompoundID]





                        print(detmethod, 'detmethod', experiment, 'experimentbla')
                        
                        if  thresholddoses[0]!=-1:
                            
                            recordingthreshold=McsPy.McsData.RawData(folder_path+'/'+thresholdfilename)
                            
                        else:
                            recordingthreshold=[]



                        D=SpikeDetection_EStim(save_to, lof, recording,
                                        CompoundID, MetadataMaskPath, tree, ExperimentID, ES_Electrodes, experiment, Analyze,
                                               detmethod, E_Stim, fs, cutoff, highpassth, order, dead_time, threshold, usewelllist, negative,
                                               positive, PipNoise, noallwells, mindur, maxdur, recordingthreshold, onlyLed, detsortfeatures, 
                                               Dosethresholddoses=thresholddoses)




                        if Analyze=='True':
                            DF= pd.DataFrame(D[0])

                            DF.to_csv(str(Path(save_to)/file+'Spikes.csv'))
                            if detsortfeatures=='True':
                                Features=D[5]

                                pd.DataFrame(Features).to_csv(str(Path(save_to)/file+'Spike_Features.csv'))
                            lof=D[2]
                            save_model(lof, '23.11.23Novelty.pkl')
                            spike_wvfs=D[6]

                            np.save(str(Path(save_to)/'spik_wvfs'), spike_wvfs)


                        ES_entry[file]=D[1]
                        LED_entry[file]=led_info(Recordingfile, MetadataFilename)
                        LabBook[file]=D[3]
                        LabBook['Version']=notebook_name()
                        LabBook['Sampling rate']=fs
                        LabBook['Threshold']=[threshold]
                        LabBook['SignNegPos']=[negative, positive]
                        LabBook['Highpass']=[highpassth]
                        LabBook['Method']=[str(detmethod)]
                        LabBook['AnanlogFhighlow']=[Ahigh, Alow]

                        stamp[file]=D[4]

                        PipArtifacts[file]=D[7]

                        print(stamp, 'ina loopstamp')

                else:

                    print('no mws for i given file in main path')









    #save_model(lof, 'Novelty.pkl') #to modify auto
    from pathlib import Path

    save_dir = Path(save_to).expanduser().resolve()

    np.save(save_dir / "ES_info.npy",   ES_entry)
    np.save(save_dir / "LED_info.npy",  LED_entry)
    np.save(save_dir / "Lab_Book.npy",  LabBook)
    np.save(save_dir / "new_stamp.npy", experiments)
    np.save(save_dir / "PipArt.npy",    PipArtifacts)
    np.save(save_dir / "stamp.npy",     stamp)

    



    return  filenames, stamp, LED_entry, LabBook


# In[ ]:






def test_main(folder, file,  MetadataMaskPath):


    """Function to detect spike and write spike detections in a table for a single recording file"""
    
    
    

    D={'Timestamp [µs]':[], 'Peak Amplitudes':[], 'Duration':[], 'Channel ID':[], 'Well ID':[], 'Well Label':[],
       'Channel Label':[], 'Experiment':[], 'Dose Label':[], 'Compound ID':[], 'ES_condition':[]}
    
    
    labs={}
    lab={}
    
    experiments={}
    
    folder_path, meta_path, saves=Fold_Folders(folder)
    
    metafname=file[:-7]+'.mws'
                    
    meta_pathf=[mp for mp in meta_path if metafname in os.listdir(mp)]




    if len(meta_pathf)>0:


        MetadataFilename=meta_pathf[0]+'/'+file[:-7]+'.mws'
        
        Recordingfilename=folder_path[0]+'/'+file
        recording=McsPy.McsData.RawData(Recordingfilename)
        #print(recording.recordings[0].event_streams.keys())

        

        tree = etree.parse(MetadataFilename)
        lab=lab_book(tree)
        WellID=np.load(MetadataMaskPath+'/'+'WellID.csv.npy')
        Welllabel=np.load(MetadataMaskPath+'/'+'WellLabel.csv.npy')
        Channellabel=np.load(MetadataMaskPath+'/'+'ChannelLabel.csv.npy')
        ExperimentID=tree.xpath('//ExperimentID/text()')[0]
        ES_Electrodes=tree.xpath('//ChannelID/text()')




        ChannelID=list(recording.recordings[0].analog_streams[0].channel_infos.keys())
        #print(list(recording.recordings[0].analog_streams[0].channel_infos.keys()))

        ChannelID=[int(i) for i in ChannelID]





        for i in ChannelID:


                index=int(recording.recordings[0].analog_streams[0].channel_infos[i].row_index)




        Compounds=tree.xpath('//CompoundID/text()')
        #print(Compounds)
        CompoundID=recording.comment

        Labels=Dilutions(tree, CompoundID)
        starts=recording.recordings[0].analog_streams[1].timestamp_index[:, 1].astype('int64')
        stops=recording.recordings[0].analog_streams[1].timestamp_index[:, 2].astype('int64')
        esim=recording.recordings[0].analog_streams[1].timestamp_index[:, 0].astype('int64')

        if ExperimentID not in experiments.keys():

            experiments[ExperimentID]={CompoundID:[]}
            experiments[ExperimentID][CompoundID].extend(esim[:].tolist())


        else:

            if CompoundID not in experiments[ExperimentID].keys():

                experiments[ExperimentID].update({CompoundID:[]})


                experiments[ExperimentID][CompoundID].extend(esim[:].tolist())

            else:

                experiments[ExperimentID][CompoundID].extend(esim[:].tolist())





        #print(starts, stops, esim, 'start, stops, esim')

        ES_dict=ES_info(recording, tree, experiments, WellID, CompoundID)

        for j in range(len(starts)):

            dose_label=Labels[j]
            dose_start=int(starts[j])
            dose_stop=int(stops[j])
            #print(dose_start, dose_stop, dose_label)

        ES_Electrodes=tree.xpath('//ChannelID/text()')
        info=led_info(Recordingfilename, MetadataFilename)

        lab['Labels']=Labels
        lab['starts']=starts
        lab['stops']=stops
        lab['LED stimulation']=info

        lab['ES_stimulation_electrodes']= ES_Electrodes
        lab['Compounds']= Compounds
        lab['Electrical stimulation']= ES_dict
        
        








    return labs, lab, recording

