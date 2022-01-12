"""
    Local Solver

    A Python program to extract amplitudes of the reflection of interest from window 
    
    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""

import numpy as np
import scipy.io as spio
from devito import *
from examples.seismic import *
from scipy.signal import hilbert
from scipy import interpolate
from calc_phase import *


def window_shot(vx_wf, vx_extrapolate, x_arr, time_index, window_length):
    """
    Windowed the shotgather

    Parameters:

    * vx_wf: 2d-arrays
        Arrays with traces, shotgather. Approximate nx ~ 50 traces, nt ~ 1700 time steps
        
    * vx_extrapolate: 2d-arrays
        Arrays with traces, shotgather. Approximate nx ~ 50 traces, nt ~ 1700 time steps

    * x_arr:1D array in x direction. 
    	Array with receiver positions inside the window  
    
    * time_index: 1D array with start time index.
 
    Returns:

    * windowed_shotgather_true_vx : 2d-arrays Arrays with traces, shotgather. Approximate nx ~ 50 traces, nt ~ 1700 time steps
    * windowed_shotgather_extr_vx : 2d-arrays Arrays with traces, shotgather. Approximate nx ~ 50 traces, nt ~ 1700 time steps

    """
    window_vx_tr1=[]
    window_vx_ex1=[]
    
    for ir in x_arr:
    	window_true_vx=vx_wf[:,ir]
    	window_vx_tr1.append(window_true_vx)
    	window_extra_vx=vx_extrapolate[:,ir]
    	window_vx_ex1.append(window_extra_vx)
    window_vx_tr1=np.asanyarray(window_vx_tr1)
    window_vx_ex1=np.asanyarray(window_vx_ex1)
    
    window_vx2=[]
    window_ex2=[]
    for i in range (len(time_index)):
    	windowed_true_vx2=window_vx_tr1.T[int(time_index[i]):int(time_index[i]+window_length),i]
    	window_vx2.append(windowed_true_vx2)
    	windowed_extr_vx2=window_vx_ex1.T[int(time_index[i]):int(time_index[i]+window_length),i]
    	window_ex2.append(windowed_extr_vx2)
    windowed_shotgather_true_vx=np.asanyarray(window_vx2)
    windowed_shotgather_extr_vx=np.asanyarray(window_ex2) 
    return   windowed_shotgather_true_vx, windowed_shotgather_extr_vx

		
def amplitude_extraction(windowed_shotgather_true_vx,windowed_shotgather_extr_vx):
    """
    Extract the amplitude from the windowed shotgather

    Parameters:

    * windowed_shotgather_true_vx: 2d-arrays
        Arrays with traces, shotgather.
        
    * windowed_shotgather_extr_vx: 2d-arrays
        Arrays with traces, shotgather.

    Returns:

    * ampl_true_norm_vx  : 1d-arrays Arrays with full amplitudes.
    * ampl_local_norm_vx : 1d-arrays Arrays with local amplitudes.

    """
    [nr, _]=windowed_shotgather_true_vx.shape
    ampl_true_vx    =np.zeros(nr)
    ampl_local_vx   =np.zeros(nr)
    for ir in range (nr):
    	ampl_true_vx[ir]    =np.max(np.abs(hilbert(windowed_shotgather_true_vx.T[:,ir])))
    	ampl_local_vx[ir]   =np.max(np.abs(hilbert(windowed_shotgather_extr_vx.T[:,ir])))
    ampl_true_norm_vx    =(ampl_true_vx - np.min(ampl_true_vx))/(np.max(ampl_true_vx)-np.min(ampl_true_vx))  
    ampl_local_norm_vx   =(ampl_local_vx - np.min(ampl_local_vx))/(np.max(ampl_local_vx)-np.min(ampl_local_vx)) 
    return  ampl_true_norm_vx,ampl_local_norm_vx	
	
def simple_plot(arr, name_output):

    """
    simple plotter

    Parameters:

    * arr: 2d-arrays.
        
    * name_output: str(out name)
    Returns:

    * figure

    """
    aspect=1.0
    perc=93
    nt,nr = arr.shape
    ratio = float(nr)/nt
    cmap = plt.get_cmap('gist_gray')
    imshow_aspect = ratio/aspect
    vmax=0.002
    vmin=-0.002
    plt.figure()
    im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect=imshow_aspect, vmax=vmax, vmin=vmin)  
    plt.savefig('fig/%s'%name_output)
	

##################################################################################################################
##Main code
##load the extrapolated shotgather
dict_extrapolate = spio.loadmat('output/shotgather_extrapolated.mat')
vx_extrapolate   = dict_extrapolate['vx_extrapolate']
vz_extrapolate   = dict_extrapolate['vz_extrapolate'] 

##load the full shot-gather
dict_extrapolate =spio.loadmat('output/shotgather_wf.mat')
vx_wf  = dict_extrapolate['vx_wf']
vz_wf  = dict_extrapolate['vz_wf']

##define time index, selection of 31 time indexes that show the interes reflection only, between trace 10 and 40.
## for the two layer model. If the model or time is modify the time index should be modify accordingly. )
time_index=[1320, 1310, 1310, 1300, 1290, 1270, 1250, 1250, 1240, 1230, 1230, 1215, 1215, 1215, 1215, 1210, 1210, 1210, 1210, 1210, 1210, 1210, 1225, 1230, 1240, 1250, 1260, 1270, 1280, 1300, 1310  ]


#############################################################
#plot time_index vs rec position
x_arr=np.arange(10,41,1)
window_length = 160
arr=vx_wf
plt.figure()
im=plt.imshow(arr, cmap='gist_gray', aspect='auto', vmax=0.02, vmin=-0.02)
plt.plot(x_arr, time_index, '*m')
plt.colorbar()	
plt.savefig('fig/QC_shot_timeindex')

#############################################################
#window the shotgather 
windowed_shotgather_true_vx, windowed_shotgather_extr_vx = 	window_shot(vx_wf, vx_extrapolate, x_arr, time_index, window_length)
#############################################################

#amplitude calculation
ampl_true_norm_vx,ampl_local_norm_vx = 	amplitude_extraction(windowed_shotgather_true_vx,windowed_shotgather_extr_vx)
#############################################################
simple_plot(arr=windowed_shotgather_true_vx.T,name_output="windowed_shotgather_true_vx.png")
simple_plot(arr=windowed_shotgather_extr_vx.T, name_output='windowed_shotgather_extr_vx.png')
###############################################
plt.figure()
im=plt.plot(ampl_true_norm_vx, '--k',label='Full', lw=2)
plt.plot(ampl_local_norm_vx, '*--r',label='local', lw=2)
plt.ylabel('Amplitude')
plt.xlabel('Receivers')
plt.legend()
plt.savefig('fig/amplitude.png')

##############################################
##Phase calculation
best_phases_true_vx, corr_coeffs_true_vx   = calc_phase(windowed_shotgather_true_vx.T, ret_corr_coeffs=True)
best_phases_local_vx,corr_coeffs_local_vx  = calc_phase(windowed_shotgather_extr_vx.T, ret_corr_coeffs=True)

##############################################
plt.figure()
im=plt.plot(best_phases_true_vx, '*--k',label='Full',   lw=2)
im=plt.plot(best_phases_local_vx,'>--r', label='Local', lw=3)
plt.title('best_phases_true_vx')
plt.ylabel('Phase (Degree)')
plt.xlabel('Receivers')
plt.legend()
plt.savefig('fig/Phase.png')
##############################################
####Save phases
outdict = dict()
outdict['best_phases_true_vx']   = best_phases_true_vx
outdict['best_phases_local_vx']  = best_phases_local_vx
outdict['corr_coeffs_true_vx']   = corr_coeffs_true_vx
outdict['corr_coeffs_local_vx']  = corr_coeffs_local_vx
spio.savemat("output/phases.mat", outdict)
##############################################

 

