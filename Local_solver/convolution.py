"""
    Local Solver

    A Python program to extrapolate the scattered wavefield to the receivers positions
    
    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""


import numpy as np
import scipy.io as spio
from devito import *
from examples.seismic import *
from scipy import signal
from scipy.signal import fftconvolve

__all__ = ['convolve_traces', 'convolve_trace_pairs', 'convolution']

def convolve_traces(trace_1, trace_2): 
    """
    Parameters:

    * Trace 1: 1D-array e.g Green function

    * Trace 2: 1D-array e.g. source wavelet 
    
    Traces don't need the same length, but I assume the 'dt' is the same

    Returns:

    * convolved trace 

    """
    return fftconvolve(trace_1, trace_2)


def convolve_trace_pairs(gather1, gather2):
    """
    Parameters:

    * gather 1: 2D-array e.g Green function

    * gather 2: 2D-array e.g. source wavelet 

    Returns:

    * convolved trace 

    """
    ntr1 = gather1.shape[1]
    ntr2 = gather2.shape[1]
    if ntr1 != ntr2:
        raise Exception("Not equal number of traces in both gathers") 

    for i in range(ntr1):
        trace1 = gather1[:,i]
        trace2 = gather2[:,i]
        outtrace = convolve_traces(trace1, trace2)
        if i == 0:
            nt_out = outtrace.size
            out    = np.zeros((nt_out, ntr1))
             
        out[:,i] = outtrace

    return out
    
def convolution(geom, wave):
    """
    Parameters:

    * geom: address to input geometry information.

    * wave: address to load scattered wavefield. 
    
    Returns:

    * Extrapolated shotgather
    * convolution output from integral  

    """
    #################################################
    #load geometry information
    ################################################
    dict_geom = spio.loadmat(geom)
    ################################################
    rec_coord    = dict_geom['rec_coord']
    src_coord    = dict_geom['src_coord']
    se_coord     = dict_geom['se_coord']
    space        = dict_geom['spacing']
    dt           = dict_geom['dt'][0][0]
    source_data  = dict_geom['source_data'][0]
    
    #####################################################
    ##load the scatter wavefield
    #####################################################
    
    dict_green_sc =spio.loadmat(wave)
    
    tzz_sc = dict_green_sc['tzz']
    txx_sc = dict_green_sc['txx']
    txz_sc = dict_green_sc['txz']
    vx_sc  = dict_green_sc['vx']
    vz_sc  = dict_green_sc['vz'] 
    #################################################
    #Define parameters
    dx  = space[0][0]
    dz  = space[0][1]
    dttx = dt
    dttz = dt
    src_coordinates=rec_coord
    #################################################
    #distance between se points
    d_se = ((se_coord[1][0]-se_coord[0][0]))/(dx*dz) 
    ##################################################
    vx_out =[]
    vz_out =[]
    for isource in src_coordinates:
    
    	ns = len(src_coordinates)
    	##load the background Green function for Vx
    	##load velocity and stress Greens tensors from the observation point (receivers positions)
    	##to the point of integration (Se) computed in unaltered model.
    	
    	dict_green_rec_se_x =spio.loadmat('output/G_Vx/Component_source_greens_se%.2f'%(isource[0]))
    	Gtzz_se_x = dict_green_rec_se_x['tzz']
    	Gtxx_se_x = dict_green_rec_se_x['txx']
    	Gtxz_se_x = dict_green_rec_se_x['txz']
    	Gvx_se_x  = dict_green_rec_se_x['vx']
    	Gvz_se_x  = dict_green_rec_se_x['vz'] 
    	
    	##Convolve traces integral vx (integral part related to x in equation (4))
    	integrand_vx  = (convolve_trace_pairs(Gvx_se_x,txz_sc)*dttx-convolve_trace_pairs(Gtxz_se_x,vx_sc)*dttx)+(convolve_trace_pairs(Gvz_se_x,tzz_sc)*dttx-convolve_trace_pairs(Gtzz_se_x,vz_sc)*dttx)
    	integral_vx   = np.sum(integrand_vx,axis=1)*d_se
    	vx_out.append(integral_vx)
    	
    	##load the background Green function for Vz
    	##load velocity and stress Greens tensors from the observation point (receivers positions)
    	##to the point of integration (Se) computed in unaltered model.
    	
    	dict_green_rec_se_z =spio.loadmat('output/G_Vz/Component_source_greens_se%.2f'%(isource[0]))
    	Gtzz_se_z = dict_green_rec_se_z['tzz']
    	Gtxx_se_z = dict_green_rec_se_z['txx']
    	Gtxz_se_z = dict_green_rec_se_z['txz']
    	Gvx_se_z  = dict_green_rec_se_z['vx']
    	Gvz_se_z  = dict_green_rec_se_z['vz']
    	
    	##Convolve traces integral vz (integral part related to z in equation (4))
    	integrand_vz = (convolve_trace_pairs(Gvx_se_z,txz_sc)*dttz-convolve_trace_pairs(Gtxz_se_z,vx_sc)*dttz)+(convolve_trace_pairs(Gvz_se_z,tzz_sc)*dttz-convolve_trace_pairs(Gtzz_se_z,vz_sc)*dttz)
    	integral_vz   = np.sum(integrand_vz,axis=1)*d_se
    	vz_out.append(integral_vz)
    	
    vx_out_f=np.asarray(vx_out)
    vz_out_f=np.asarray(vz_out)
    print('integral done')
    ################################
    
    # load velocity background.
    dict_green_back =spio.loadmat('output/shotgather_src_nonp_%.2f'%(src_coord[0][0]))
    vx_rec_back  = dict_green_back['vx']
    vz_rec_back  = dict_green_back['vz'] 
    #######################################
    
    #convolve background with source 
    nt,nr    = vz_rec_back.shape
    print('nt',nt, 'nr', nr)
    
    vx_out_b = np.zeros((nt, nr))
    for ir in range (vx_rec_back.shape[1]):
    	vx_rec_src = convolve_traces(vx_rec_back[:,ir],source_data)
    	vx_out_b[:,ir]  = vx_rec_src[:nt]*(dt) 
    	
    vz_out_b = np.zeros((nt, nr))
    for ir in range (vz_rec_back.shape[1]):
    	vz_rec_src = convolve_traces(vz_rec_back[:,ir],source_data)
    	vz_out_b[:,ir]  = vz_rec_src[:nt]*(dt)  
    	
    # add background to the integral (full equation 4)
    vx_extrapolate = vx_out_b[:]+ vx_out_f[:, 2:nt+2].T
    vz_extrapolate = vz_out_b[:]+ vz_out_f[:, 2:nt+2].T
    
    ##########################################################################################
    #save outputs
    #############################################################################################
    
    outdict = dict()
    outdict['vx_out_integral']= vx_out_f
    outdict['vz_out_integral']= vz_out_f
    outdict['vx_extrapolate']= vx_extrapolate
    outdict['vz_extrapolate']= vz_extrapolate
    spio.savemat('output/shotgather_extrapolated.mat', outdict)
    
    print('extrapolation done')
    print('convolution done')
    extrapolate_shotgather = {'vx_extrapolate':vx_extrapolate, 'vz_extrapolate':vz_extrapolate}
    return extrapolate_shotgather
