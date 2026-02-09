#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 17:12:14 2024

@author: jaramillo
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy.signal as sig
import scipy.integrate as integrate
import scipy.stats as sc_stats
import seaborn as sns
import os
import get_CCG_functions as gcg
##### IMPORT files from jupyter notebook/NWB

os.chdir('/Users/jaramillo/map-ephys/notebook/workshop/April2024NWB')
def load_onesession(sub_id, session, regions): 
    """
    

    Parameters
    ----------
    sub_index : index for subject
    session_index : index for session
    Returns
    -------
    allspikes, allunits, stats
    
    """
    path = '/Users/jaramillo/map-ephys/notebook/workshop/April2024NWB/NWBdata'
    os.chdir(path+'/'+'sub-'+sub_id+'/analysis')    
    file0 = regions[0] +'_notrials'+'sub-'+str(sub_id)+'_ses-'+str(session)+'_allunits.pkl'
    file1 = 'locationsCCF'+ 'sub-'+str(sub_id)+'_ses-'+str(session)+'_allunits.pkl'
  
    with open(file0, 'rb') as f:  # open a text file
        allspikes = pickle.load(f) # # 
    with open(file1, 'rb') as f:  # open a text file
        allCCF = pickle.load(f) # # 
    with open('CCG'+regions[0]+'-'+regions[1]+'sub-'+str(sub_id)+'_ses-'+str(session)+'.pkl', 'rb') as f:  # open a text file
        svcorr_vec = pickle.load(f)
    
    return allspikes, allCCF, svcorr_vec


def get_allCCG(regions, alldirectories):#p
    '''
    
    regions is a list of two strings that define the cross correlation e.g., ['left ALM', 'left ALM']
            '''
    path = '/Users/jaramillo/map-ephys/notebook/workshop/April2024NWB/NWBdata'
    os.chdir(path)
    if alldirectories =='yes':
       directories = os.listdir
    else:
        directories = ['sub-484674', 'sub-480133', 'sub-480135', 'sub-480928', 'sub-484672', 'sub-440959', 'sub-480927',  'sub-484675', 'sub-484673', 'sub-480134']#['sub-455220']
    for dir1 in directories:
        if not dir1.startswith('.'):
            sessions = gcg.get_sessions(path, dir1)
            sub_id = dir1[4:]
            for session in sessions: ### loops over sessions within a subdirectory 
                os.chdir(path+'/'+dir1+'/analysis')
                spikes_pooled = {}
                file0 = regions[0] +'_notrials'+'sub-'+str(sub_id)+'_'+str(session)+'_allunits.pkl'
                file1 = regions[1] +'_notrials'+'sub-'+str(sub_id)+'_'+str(session)+'_allunits.pkl'
                if os.path.isfile(file0)==True and os.path.isfile(file1)==True :
                    print('session '+ (sub_id+session)+' crosscorrelation between ' +regions[0]+ ' and '+ regions[1])
                    with open(file0, 'rb') as f:  # open a text file
                        spikes_pooled[regions[0]] = pickle.load(f) # # 
                    with open(file1, 'rb') as f:  # open a text file
                        spikes_pooled[regions[1]] = pickle.load(f) # # 
                    sparse1, sparse2 = gcg.getsparsematrix(spikes_pooled[regions[0]], spikes_pooled[regions[1]])
                    corr_vec, filt_time, ALM_FR, Thal_FR = gcg.cross_corr_sam(sparse1, sparse2)
                    with open('CCG'+regions[0]+'-'+regions[1]+'sub-'+str(sub_id)+'_'+str(session)+'.pkl', 'wb') as f:  # open a text file
                        pickle.dump(corr_vec, f) # 
                    del corr_vec
                else:
                    print('session '+ (sub_id+session)+' does not have'+regions[0]+ ' and '+ regions[1])
                    
def get_allselectivityvectors(regions, *params):
     ###load parameters
     peak_th, norm, peakstrength, strict_contraipsi, alldirectories = params
     ### crosscorrelation parameters
    
     ### variable initialization
     
     all_sel_cx = np.zeros((1,2))
     all_sel_th = np.zeros((1,2))
     all_CCF = []#np.zeros((1,1,2))
     all_CCF_cx = np.zeros((1,2))
     all_CCF_th = np.zeros((1,2))

    
     ####spike vec parameters
     binsize = 0.001#dt
     timevec = np.arange(-2.5,1.5,binsize)
     #### load data
     path = '/Users/jaramillo/map-ephys/notebook/workshop/April2024NWB/NWBdata'
     os.chdir(path)
     #alldirectories = 'no'
     if alldirectories =='yes':
        directories = os.listdir()
     else:
         directories = ['sub-456772', 'sub-484673']#,'sub-455219', 'sub-480927']
     for dir1 in directories:
         if not dir1.startswith('.'):
             sessions = gcg.get_sessions(path, dir1)
             sub_id = dir1[4:]
             for session in sessions: ### loops over sessions within a subdirectory 
                 os.chdir(path+'/'+dir1+'/analysis')
                 #print('dir1', dir1)
                 spikes_seg = {}
                 CCF_allregions =  {}
                 file0 = regions[0] +'_withtrials'+'sub-'+str(sub_id)+'_'+str(session)+'_allunits.pkl'
                 file1 = regions[1] +'_withtrials'+'sub-'+str(sub_id)+'_'+str(session)+'_allunits.pkl'
                 hemi = regions[0][0:4]
                 #print('file0', file0)
                 #print('file1', file1)

                 if os.path.isfile(file0)==True and os.path.isfile(file1)==True :
                     with open(file0, 'rb') as f:  # open a text file
                         spikes_seg[regions[0]] = pickle.load(f) # # 
                     with open(file1, 'rb') as f:  # open a text file
                         spikes_seg[regions[1]] = pickle.load(f) # # 
                     with open('session'+str(session)[4:]+'_sub'+str(sub_id)+'_stats.pkl', 'rb') as f: 
                         stats= pickle.load(f) 
                     with open('locationsCCF'+'sub-'+regions[0]+str(sub_id)+'_'+str(session)+'_allunits.pkl', 'rb') as f:  
                         CCF_allregions[regions[0]] = pickle.load(f) # # 
                     with open('locationsCCF'+'sub-'+regions[1]+str(sub_id)+'_'+str(session)+'_allunits.pkl', 'rb') as f:  
                         CCF_allregions[regions[1]] = pickle.load(f) # # 
                         
                         
                     print(str(sub_id)+str(session))
                     spikevec_ALM, spikevec_Thal, sel_vec_cx, sel_vec_th, ALM_FR, Thal_FR  = gcg.spikestosel(spikes_seg[regions[0]], spikes_seg[regions[1]], timevec, stats, hemi, 0.05)
                     #### information about number of selective neurons for a given session, and cummulatively
                     if not sel_vec_cx:
                         print(str(sub_id)+str(session) + 'has no selvec')

                     sel_vec_cx_array = np.array(sel_vec_cx)
                     sel_vec_th_array = np.array(sel_vec_th)
                     #print((sel_vec_cx_array))
                     all_sel_cx = np.vstack([all_sel_cx, sel_vec_cx_array])
                     all_sel_th = np.vstack([all_sel_th, sel_vec_th_array])
                     all_CCF.append(CCF_allregions)
     for j in range(0, len(all_CCF)):
           all_CCF_cx = np.vstack([all_CCF_cx, all_CCF[j][regions[0]]])   
           all_CCF_th = np.vstack([all_CCF_th, all_CCF[j][regions[1]]])   

     return all_sel_cx, all_sel_th, all_CCF_cx, all_CCF_th, all_CCF            
                     
                            

def get_allpeaks(regions,all_sessiondic,*params):
    '''
    Calculates the peaks of the crosscorrelation function

    Parameters
    ----------
    regions : list of strings, defines the regions to import crosscorrelation data from
        
    *params : list.     peak_th, norm, peakstrength, strict_contraipsi = params


    Returns
    -------
    all_contrapeaks : accumulates all contrapeaks as defined by the params
    all_ipsipeaks : accumulates all ipsipeaks as defined by the params
    all_nonselpeaks : ccumulates all nonselective peaks as defined by the params

    '''
    ###load parameters
    peak_th, norm, peakstrength, strict_contraipsi, alldirectories = params
    ### crosscorrelation parameters
    dt = 0.0001
    maxlag = 20e-3
    Nlag = int(maxlag/dt)
    filt_time = dt*np.arange(-Nlag, Nlag+1)
    ### variable initialization
    all_contrapeaks = []
    all_ipsipeaks = []
    all_nonselpeaks = []
    all_sel_pre = []
    all_sel_post = []
    all_PC_pre = []
    all_PC_post = []
    all_efficacy = []
    all_CCG = np.zeros((0, len(filt_time)))
    all_peakindices = np.zeros(2)
    all_CCF = np.zeros((1,2,3))
    all_sessions = [] 

    totalcontra_cx = 0
    totalipsi_cx = 0
    totalcontra_th = 0
    totalipsi_th = 0
    ####spike vec parameters
    binsize = 0.001#dt
    timevec = np.arange(-2.5,1.5,binsize)
    #### load data
    path = '/Users/jaramillo/map-ephys/notebook/workshop/April2024NWB/NWBdata'
    os.chdir(path)
    #alldirectories = 'no'
    if alldirectories =='yes':
       directories = os.listdir()
    else:
        directories = ['sub-455219', 'sub-456772']
    for dir1 in directories:
        if not dir1.startswith('.'):
            sessions = gcg.get_sessions(path, dir1)
            sub_id = dir1[4:]
            for session in sessions: ### loops over sessions within a subdirectory 
                os.chdir(path+'/'+dir1+'/analysis')
                print('dir1', dir1)
                spikes_seg = {}
                CCF_allregions =  {}
                file0 = regions[0] +'_withtrials'+'sub-'+str(sub_id)+'_'+str(session)+'_allunits.pkl'
                file1 = regions[1] +'_withtrials'+'sub-'+str(sub_id)+'_'+str(session)+'_allunits.pkl'
                hemi = regions[0][0:4]
                #print('file0', file0)
                #print('file1', file1)

                if os.path.isfile(file0)==True and os.path.isfile(file1)==True :
                    #print('session '+ (sub_id+session)+' crosscorrelation between ' +regions[0]+ ' and '+ regions[1])
                    with open(file0, 'rb') as f:  # open a text file
                        spikes_seg[regions[0]] = pickle.load(f) # # 
                    with open(file1, 'rb') as f:  # open a text file
                        spikes_seg[regions[1]] = pickle.load(f) # # 
                    with open('session'+str(session)[4:]+'_sub'+str(sub_id)+'_stats.pkl', 'rb') as f: 
                        stats= pickle.load(f) 
                    with open('CCG'+regions[0]+'-'+regions[1]+'sub-'+str(sub_id)+'_'+str(session)+'.pkl', 'rb') as f:  # open a text file
                        corr_vec = pickle.load(f)
                    with open('locationsCCF'+'sub-'+regions[0]+str(sub_id)+'_'+str(session)+'_allunits.pkl', 'rb') as f:  
                        CCF_allregions[regions[0]] = pickle.load(f) # # 
                    with open('locationsCCF'+'sub-'+regions[1]+str(sub_id)+'_'+str(session)+'_allunits.pkl', 'rb') as f:  
                        CCF_allregions[regions[1]] = pickle.load(f) # # 
                        
                        
                    #print("load CCG sub"+str(sub_id)+" ses "+str(session)+"complete")
                    
                    #spikevec_ALM, spikevec_Thal, sel_vec_cx, sel_vec_th, ALM_FR, Thal_FR  = gcg.spikestosel(spikes_seg[regions[0]], spikes_seg[regions[1]], timevec, stats, hemi, 0.05)
                    sel_vec_cx = all_sessiondic[file0+' _sel'] 
                    sel_vec_th  = all_sessiondic[file1+' _sel'] 
                    ALM_FR = all_sessiondic[file0+' _FR']  
                    Thal_FR = all_sessiondic[file1+' _FR'] 
                    #### information about number of selective neurons for a given session, and cummulatively
                    """
                    sel_vec_cx_array = np.array(sel_vec_cx)
                    sel_vec_th_array = np.array(sel_vec_th)
                  
                    totalcontra_cx+=np.sum(sel_vec_cx_array[:,0]=='contra')
                    totalipsi_cx+=np.sum(sel_vec_cx_array[:,0]=='ipsi')
                    totalcontra_th+=np.sum(sel_vec_th_array[:,0]=='contra')
                    totalipsi_th+=np.sum(sel_vec_th_array[:,0]=='ipsi')
                    """
                    
                    ccgs_norm = gcg.CCG_norm(corr_vec, filt_time, ALM_FR, Thal_FR, norm)
                    session_sub = session+str(sub_id)
                    #print('session_sub', session_sub)
                    #print('session',session)
                    peakindices_alltrials, CCF_alltrials, ccgwithpeak_alltrials, allpeaks_alltrials, peak_sel_alltrials, allcounters_alltrials, sessions_alltrials \
                        = gcg.peak_filt(dt, filt_time, ccgs_norm, sel_vec_cx, sel_vec_th, CCF_allregions[regions[0]], CCF_allregions[regions[1]],ALM_FR, Thal_FR,peak_th, session_sub, 0.015, peakstrength)        
                    
                    del corr_vec
                    print('sel_vec_cxshape', np.shape(sel_vec_cx))
                    print('peakshape', np.shape(peak_sel_alltrials))
                    print('CCGshape', np.shape(ccgwithpeak_alltrials))

                    if len(peak_sel_alltrials)>0:
                        sel_pre, sel_post, efficacy = np.transpose(peak_sel_alltrials)

                    else:
                        sel_pre, sel_post, efficacy = 'n', 'n','n'
                    
                    """
                    
                    contra_peaks = allpeaks_alltrials['peaks_contracontra'] + strict_contraipsi*allpeaks_alltrials['peaks_contranon']
                    ipsi_peaks = allpeaks_alltrials['peaks_ipsiipsi'] + strict_contraipsi*allpeaks_alltrials['peaks_ipsinon']
                    nonsel_peaks = allpeaks_alltrials['peaks_nonnon'] + allpeaks_alltrials['peaks_ipsinon'] + allpeaks_alltrials['peaks_contranon'] + allpeaks_alltrials['peaks_mixed']
                    all_contrapeaks = np.hstack([all_contrapeaks, contra_peaks])
                    all_ipsipeaks = np.hstack([all_ipsipeaks, ipsi_peaks])
                    all_nonselpeaks = np.hstack([all_nonselpeaks, nonsel_peaks])
                    """
                    if len(peak_sel_alltrials)>0:
                        all_CCG = np.vstack([all_CCG, ccgwithpeak_alltrials])
                        all_peakindices = np.vstack([all_peakindices, peakindices_alltrials])#all_peakindices.append(peakindices_alltrials)
                        print('shape_allpeakindices', np.shape(all_peakindices))
                        #print('CCFshape', np.shape(CCF_alltrials))
                        #print('CCF', (CCF_alltrials))
                        
                        CCF_alltrials = np.array(CCF_alltrials)
                        all_CCF = np.vstack([all_CCF, CCF_alltrials])
                        all_sessions.append(sessions_alltrials)# = np.vstack([all_sessions, sessions_alltrials])
                        all_sel_pre = np.hstack([all_sel_pre, sel_pre])
                        all_sel_post = np.hstack([all_sel_post, sel_post])
                        all_efficacy= np.hstack([all_efficacy, efficacy])

    all_CCF = all_CCF[1:,:,:]
    all_peakindices = all_peakindices[1:, :]   
    all_sessions = [x for xs in all_sessions for x in xs]
               
    totalsel  = [totalcontra_cx, totalipsi_cx, totalcontra_th, totalipsi_th]
    """plt.figure()
    sns.kdeplot(all_contrapeaks, fill = 'True', color = 'steelblue').set(xlim=0)     
    sns.kdeplot(all_ipsipeaks, fill = 'True', color = 'red').set(xlim=0)
    plt.figure()
    #plt.hist(all_contrapeaks, bins=np.arange(min(all_contrapeaks), max(all_contrapeaks) + binwidth, binwidth), color = 'steelblue')    
    #plt.hist(all_ipsipeaks, bins=np.arange(min(all_ipsipeaks), max(all_ipsipeaks) + binwidth, binwidth), color='red')
    #print('mean_contra', np.mean(all_contrapeaks))    
    #print('mean_ipsi', np.mean(all_ipsipeaks))
    plt.figure()
    #plt.scatter(all_sel_pre, all_sel_post, s = 20000*all_efficacy)
    #plt.axvline(x = 0, color = 'k')
    #plt.axhline(y = 0, color = 'k')
    #plt.axvline(x = 0, color = 'k')
    #plt.axhline(y = 0, color = 'k')
    """
    return all_sel_pre, all_sel_post, all_efficacy, all_CCG, totalsel, all_peakindices, all_CCF, all_sessions, filt_time         

#all_contrapeaks, all_ipsipeaks, all_nonselpeaks, all_sel_pre, all_sel_post, all_efficacy, all_CCG, totalsel, all_peakindices, all_CCF, all_sessions, filt_time         

def plot_EAB(new, all_sessiondic_L, all_sessiondic_R,  *params):
    
    if new==1:
       params_CCG = params
       all_sel_pre_cxth_L, all_sel_post_cxth_L, all_efficacy_cxth_L, ccgwithpeak_alltrials_cxth_L, totalselL , all_peakindices_L, all_CCF_bad_L, all_sessions_L, filt_time = get_allpeaks(['left ALM', 'left Thalamus'], all_sessiondic_L,*params_CCG)  
       all_sel_pre_cxth_R, all_sel_post_cxth_R, all_efficacy_cxth_R, ccgwithpeak_alltrials_cxth_R, totalselR , all_peakindices_R, all_CCF_bad_R, all_sessions_R, filt_time = get_allpeaks(['right ALM', 'right Thalamus'], all_sessiondic_R,*params_CCG)        

    else: 
        
       with open('dataEABLoct21'+'.pkl', 'rb') as f:  # open a text file
           all_contrapeaks_cxth_L, all_ipsipeaks_cxth_L, all_nonselpeaks_cxth_L, all_sel_pre_cxth_L, all_sel_post_cxth_L, all_efficacy_cxth_L, ccgwithpeak_alltrials_cxth_L, totalselL , all_peakindices_L, all_CCF_bad_L, filt_time \
         = pickle.load(f)
       with open('dataEABRoct21'+'.pkl', 'rb') as f:  # open a text file
            all_contrapeaks_cxth_R, all_ipsipeaks_cxth_R, all_nonselpeaks_cxth_R, all_sel_pre_cxth_R, all_sel_post_cxth_R, all_efficacy_cxth_R, ccgwithpeak_alltrials_cxth_R, totalselR , all_peakindices_R, all_CCF_bad_R, filt_time \
          = pickle.load(f)
        
    allunitindices_L= [[503,649], [502,649], [193,647], [163,649], [163,648], [153,649], [152,1009], [152,653], [152,649], [152,648], [146,1009], [146,649], [145,649], [128,649], [502, 1604], [390, 1846], [256, 605], [251, 605], [57,1480], [331, 599], [256, 546], [365, 1722], [320, 1250], [390,1634], [220,466], [256,543]]
    allunitindices_R = [[1190,2016], [1020,2293], [1020, 2073], [1010, 2073], [998, 2293], [998,2097], [998, 2091], [998, 2089], [989, 2097], [361, 1009], [285, 536], [236, 588], [152, 719], [54, 718], [1120, 2577], [1120, 2265], [730, 2649]]
    index_L = gcg.get_indices_manualCCG(all_peakindices_L, allunitindices_L)
    index_R = gcg.get_indices_manualCCG(all_peakindices_R, allunitindices_R)
    all_sel_pre_cxth = np.hstack([all_sel_pre_cxth_L[index_L], all_sel_pre_cxth_R[index_R]])
    all_sel_post_cxth = np.hstack([all_sel_post_cxth_L[index_L], all_sel_post_cxth_R[index_R]])
    all_efficacy_cxth= np.hstack([all_efficacy_cxth_L[index_L],all_efficacy_cxth_R[index_R]])
    gcg.plots_histogramCCG(all_sel_pre_cxth,all_sel_post_cxth, all_sel_pre_cxth, all_sel_post_cxth, all_efficacy_cxth, radius = 2000.0)               
 
params_CCG = [1,"pre", "integral", 0, 'yes' ]#peak_th, norm, peak_strength, strict_contraipsi, alldirectories = params


#get_allCCG(['left ALM', 'left Thalamus'], 'no')
#get_allCCG(['left Thalamus', 'right Thalamus'], 'no')
#get_allCCG(['right ALM', 'right Thalamus'], 'no')


#all_contrapeaks_cx_L, all_ipsipeaks_cx_L, all_nonselpeaks_cx_L, all_sel_pre_cx_L, all_sel_post_cx_L, all_efficacy_cx_L, ccgwithpeak_alltrials_cx_L, totalselLcx , all_peakindices_Lcx, all_CCF_bad_L, filt_time  = get_allpeaks(['left ALM', 'left ALM'], *params_CCG)
#all_contrapeaks_cxth_R, all_ipsipeaks_cxth_R, all_nonselpeaks_cxth_R, all_sel_pre_cxth_R, all_sel_post_cxth_R, all_efficacy_cxth_R, ccgwithpeak_alltrials_cxth_R, totalselR , all_peakindices_R, all_CCF_bad_R, filt_time = get_allpeaks(['right ALM', 'right Thalamus'], *params_CCG)
#all_contrapeaks_cxth_L, all_ipsipeaks_cxth_L, all_nonselpeaks_cxth_L, all_sel_pre_cxth_L, all_sel_post_cxth_L, all_efficacy_cxth_L, ccgwithpeak_alltrials_cxth_L, totalselL , all_peakindices_L, all_CCF_bad_L, filt_time = get_allpeaks(['left ALM', 'left Thalamus'], *params_CCG)




###### allspikes contains all spike times for all units. allunits is segmented into trials, as defined by the export file
    
    
### criteria/parameters
## before export: hit vs miss, all trial types? 
## after export
### selectivity: epoch where selectivity is calculated, p value
### CCG calculation: timebin, lag. 
### peakindices and CCG analysis:    peak to std, height vs integral, normalization (e.g., pre vs post)  
### definition of contra/ipsi peak? 
    
    
#### obtain selectivity and firing rate from spike trains segmented into trials
#timevec, spikevec_ALM, spikevec_Thal, sel_vec_cx, sel_vec_th, ALM_FR, Thal_FR  = spikestosel(ALM_units, Thal_units, 'all', stats, 'left', 0.05)

#### obtain spiketrain matrix (time, units), from spike trains that are not segmented intro trials
#sparse, sparse_th = getsparsematrix(allspikes_ALM, allspikes_Thal)

### obtain cross correlogram from sparse matrices for cortex and/or thalamus. ccgs_norm applies a normalization, pre or postsynaptic
#ccgs_allunits, filt_time, ALM_FR, Thal_FR = cross_corr_sam(sparse, sparse)  
#ccgs_norm = CCG_norm(ccgs, filt_time, ALM_FR, Thal_FR, 's')

### obtain different peaks
#peakindices_alltrials, ccgwithpeak_alltrials, allpeaks_alltrials, peak_sel_alltrials, allcounters_alltrials = peak_filt(0.0001, filt_time, ccgs_norm_allunits, sel_vec_cx, sel_vec_cx,ALM_FR, ALM_FR,5.5, 0.015, 'cxandcx')#allpeaksipsi = allpeaks_alltrials['peaks_ipsiipsi'] + allpeaks_alltrials['peaks_ipsinon']
#allpeakscontra = allpeaks_alltrials['peaks_contracontra'] + allpeaks_alltrials['peaks_contranon']


### plots
