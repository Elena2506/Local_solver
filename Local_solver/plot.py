"""
    Local Solver

    A Python program to plot 
    
    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""

import numpy as np
import scipy.io as spio 
from devito import *
from examples.seismic import *
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from convolution import convolve_traces
__all__ = ['plot_injection', 'plot_model', 'plot_shotgather','plot_traces', 'simple_model_plot', 'plot_simple_shotgaher', 'plot_simple_traces']



def plot_injection(geom, v_local, model, mo):
	#################################################################################
	####### Load geometry full
	#################################################################################
	dict_geom =spio.loadmat(geom)
	v1_f      = dict_geom['v1']
	src_coord = dict_geom['src_coord']
	nbl       = dict_geom['nbl'][0][0]
	i=src_coord[0][0]
	
	
	################################################################################
	####### Load v local
	#################################################################################
	dict_local =spio.loadmat(v_local)
	v1_l = dict_local['v1']
	
	#################################################################################
	####### Load points defining domains, Se and S
	#################################################################################
	
	dict_fullmodel =spio.loadmat(model)
	space   = dict_fullmodel['spacing']
	xmin=dict_fullmodel['xmin'][0][0]
	xmax=dict_fullmodel['xmax'][0][0]
	zmin=dict_fullmodel['zmin'][0][0]
	zmax=dict_fullmodel['zmax'][0][0]

	xmin_inj=dict_fullmodel['xmin_inj'][0][0]
	xmax_inj=dict_fullmodel['xmax_inj'][0][0]
	zmin_inj=dict_fullmodel['zmin_inj'][0][0]
	zmax_inj=dict_fullmodel['zmax_inj'][0][0]

	x_max_local=dict_fullmodel['xmax_local'][0][0]
	x_min_local=dict_fullmodel['xmin_local'][0][0]
	z_max_local=dict_fullmodel['zmax_local'][0][0]
	z_min_local=dict_fullmodel['zmin_local'][0][0]

	xmax_inj_local=dict_fullmodel['xmax_inj_local'][0][0]
	xmin_inj_local=dict_fullmodel['xmin_inj_local'][0][0]
	zmin_inj_local=dict_fullmodel['zmin_inj_local'][0][0]
	zmax_inj_local=dict_fullmodel['zmax_inj_local'][0][0]

	###############################################################################
	### Define for full domain 
	ix, iz = v1_f.shape
	x_max=int(ix-(nbl))
	z_max=int(iz-(nbl))
	x_extent_max=int(ix-(2*nbl))*space[0][0]
	z_extent_max=int(iz-(2*nbl))*space[0][1]

	###############################################################################
	### Define for local domain
	ixl,izl=v1_l.shape
	x_max_l=int(ixl-(nbl))
	z_max_l=int(izl-(nbl))
	x_extent_max_l=int(ixl-(2*nbl))*space[0][0]
	z_extent_max_l=int(izl-(2*nbl))*space[0][1]

	itime = 970
	plt.figure(1)
	plt.subplot(211)
	##Full
	plt.imshow(v1_f[nbl:x_max,nbl:z_max].T, cmap='seismic',extent=[0,x_extent_max, z_extent_max, 0], vmax=np.max(v1_f),vmin=np.min(v1_f))
	x1l=(xmin)*space[0][0]
	x2l=(xmax)*space[0][0]
	z1l=(zmin)*space[0][1]
	z2l=(zmax)*space[0][1]

	x1l_i=(xmin_inj)*space[0][0]
	x2l_i=(xmax_inj)*space[0][0]
	z1l_i=(zmin_inj)*space[0][1]
	z2l_i=(zmax_inj)*space[0][1]
	plt.plot((x1l_i,x1l_i,x2l_i,x2l_i,x1l_i),(z1l_i,z2l_i,z2l_i,z1l_i,z1l_i),'-b',lw=5)

	plt.plot(src_coord[0][0], src_coord[0][1], 'ok')
	plt.title('Full')
	plt.xlabel('Distance x (m)')
	plt.ylabel('Depth (m)')

	##Local
	plt.subplot(212)
	plt.imshow(v1_l[nbl:x_max_l,nbl:z_max_l].T, cmap='seismic',extent=[x1l,x2l, z2l, z1l], vmax=np.max(v1_f),vmin=np.min(v1_f))
	x1_f = 0
	x2_f = x_extent_max
	z1_f = 0
	z2_f = z_extent_max

	plt.plot((x1l_i,x1l_i,x2l_i,x2l_i,x1l_i),(z1l_i,z2l_i,z2l_i,z1l_i,z1l_i),'-b',lw=5)
	plt.ylabel('Depth (m)')
	plt.xlabel('Distance x (m)')
	plt.title('Local')
	plt.tight_layout()
	plt.savefig('fig/Injection_%s_src_%.2f.png'%(mo,i))
	
	print('Finished plot injection')
	

def plot_model(model, geom):

	###load the model 
	dict_submodel =spio.loadmat(model)

	vp = dict_submodel['vp']
	vs = dict_submodel['vs']
	rho= dict_submodel['rho']
	space = dict_submodel['spacing']
	shaping = dict_submodel['shape'] 
	spacing = (space[0][0], space[0][1])

	### load index of the local domain, injection and extrapolation surface referencing the full domain  ####
	xmin=dict_submodel['xmin'][0][0]
	xmax=dict_submodel['xmax'][0][0]
	zmin=dict_submodel['zmin'][0][0]
	zmax=dict_submodel['zmax'][0][0]

	xmin_inj= dict_submodel['xmin_inj'][0][0]
	xmax_inj= dict_submodel['xmax_inj'][0][0]
	zmin_inj= dict_submodel['zmin_inj'][0][0]
	zmax_inj= dict_submodel['zmax_inj'][0][0]

	xmin_se = dict_submodel['xmin_se'][0][0]
	xmax_se = dict_submodel['xmax_se'][0][0]
	zmin_se = dict_submodel['zmin_se'][0][0]

	#######load source coordinates 
	dict_src =spio.loadmat(geom)
	src_coord=dict_src['src_coord'][0]

	##############################################
	#Define geometry parameters
	nx, nz = vp.shape
	x_min = 0.0
	z_min = 0.0
	x_max = (nx-1)*space[0][0]
	z_max = (nz-1)*space[0][1]

	x_min_km = x_min/1000.0; x_max_km = x_max/1000.0; z_min_km = z_min/1000.0; z_max_km = z_max/1000.0; 

	arr=vp.T
	plt.figure()
	im = plt.imshow(arr, extent=[x_min_km, x_max_km,z_max_km, z_min_km])
	x1l=(xmin)*spacing[0]/1000
	x2l=(xmax)*spacing[0]/1000
	z1l=(zmin)*spacing[1]/1000
	z2l=(zmax)*spacing[1]/1000
	plt.plot((x1l,x1l,x2l,x2l,x1l),(z1l,z2l,z2l,z1l,z1l),'-r',label='local domain' ,lw=2)
	x1se=(xmin_se)*spacing[0]/1000
	x2se=(xmax_se)*spacing[0]/1000
	z1se=(zmin_se)*spacing[1]/1000
	plt.plot((x1se,x2se),(z1se,z1se),'-k',label='extrapolation surface Se' ,lw=2)
	x1l=(xmin_inj)*spacing[0]/1000
	x2l=(xmax_inj)*spacing[0]/1000
	z1l=(zmin_inj)*spacing[1]/1000
	z2l=(zmax_inj)*spacing[1]/1000
	plt.plot((x1l,x1l,x2l,x2l,x1l),(z1l,z2l,z2l,z1l,z1l),'-b',label='injection' ,lw=2)
	plt.plot(src_coord[0]/1000, src_coord[1]/1000, '*r', lw=10, label='source')

	plt.colorbar(im,label='V(km/s)', fraction=0.02, pad=0.04)
	plt.title('P-velocity full' )
	plt.xlabel('Horizontal coordinate (km)')
	plt.ylabel('Depth (km)')  
	plt.savefig('fig/model-vp.png')
	print('Finished plot model')
	
		

def plot_shotgather(shotgather_extra, shotgather_pert, geom ):
	###load extrapolate wavefield
	dict_extrapolate =spio.loadmat(shotgather_extra)
	vx_out_integral  = dict_extrapolate['vx_out_integral']
	vz_out_integral  = dict_extrapolate['vz_out_integral'] 
	vx_extrapolate   = dict_extrapolate['vx_extrapolate']
	vz_extrapolate   = dict_extrapolate['vz_extrapolate']
	
	##load velocity from full model pert.
	dict_green_f =spio.loadmat(shotgather_pert)
	vx_rec_f  = dict_green_f['vx']
	vz_rec_f  = dict_green_f['vz']
	nt, nr = vx_rec_f.shape
	
	dict_geom = spio.loadmat(geom)
	source_data  = dict_geom['source_data'][0]
	dt           = dict_geom['dt'][0][0]  
	spacing      = dict_geom['spacing']
	##############################################################
	# for comparison propose convolve full with source 
	##############################################################
	vx_out = np.zeros((nt, nr))
	for ir in range (vx_rec_f.shape[1]):
		vx_rec_src = convolve_traces(vx_rec_f[:,ir],source_data) 
		vx_out[:,ir]  = vx_rec_src[:nt]*(dt) 
	vz_out = np.zeros((nt, nr))
	for ir in range (vz_rec_f.shape[1]):
		vz_rec_src = convolve_traces(vz_rec_f[:,ir],source_data) 
		vz_out[:,ir]  = vz_rec_src[:nt]*(dt)
	
	outdict = dict()
	outdict['vx_wf']= vx_out
	outdict['vz_wf']= vz_out
	spio.savemat('output/shotgather_wf.mat', outdict)
	
	#########################################
	##Plot for Vx
	#########################################
	cmap = plt.get_cmap('gist_gray') 
	vmax=np.max(vx_extrapolate)/200
	vmin=np.min(vx_extrapolate)/200
	aspect=1.0
	perc=93
	x_extent_max=(vx_extrapolate.shape[1]*10*spacing[0][0])/1000
	z_extent_max=(vx_extrapolate.shape[0]*dt)/1000
	arr = vx_out
	plt.figure()
	plt.subplot(211)
	im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect='auto', vmax=vmax, vmin=vmin, extent=[0,x_extent_max, z_extent_max, 0])
	plt.title('Full')
	plt.xlabel('Horizontal coordinate (km)')
	plt.ylabel('Time (s)')

	plt.subplot(212)
	arr = vx_extrapolate[:] 
	im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect='auto', vmax=vmax, vmin=vmin, extent=[0,x_extent_max, z_extent_max, 0])
	plt.title('Local')
	plt.xlabel('Horizontal coordinate (km)')
	plt.ylabel('Time (s)')
	plt.tight_layout()
	plt.savefig('fig/shotgathers_vx.png')
	
	#########################################
	##Plot for Vz  
	#########################################
	
	vmax=np.max(vz_extrapolate)/1000
	vmin=np.min(vz_extrapolate)/1000
	
	arr = vz_out
	plt.figure()
	plt.subplot(211)
	im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect='auto', vmax=vmax, vmin=vmin, extent=[0,x_extent_max, z_extent_max, 0])
	plt.title('Full')
	plt.xlabel('Horizontal coordinate (km)')
	plt.ylabel('Time (s)')

	plt.subplot(212)
	arr = vz_extrapolate[:] 
	im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect='auto', vmax=vmax, vmin=vmin, extent=[0,x_extent_max, z_extent_max, 0])
	plt.title('Local')
	plt.xlabel('Horizontal coordinate (km)')
	plt.ylabel('Time (s)')
	plt.tight_layout()
	plt.savefig('fig/shotgathers_vz.png')
	
	#########################################
	##vx, vz out of convolution
	
	aspect=1.0
	perc=93
	arr = vx_out_integral[:].T 
	nt,nr = arr.shape
	ratio = float(nr)/nt
	imshow_aspect = ratio/aspect

	plt.figure()
	im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect=imshow_aspect, vmax=vmax, vmin=vmin)
	plt.colorbar()
	plt.title('vx_conv')
	plt.xlabel('nx')
	plt.ylabel('nt')
	plt.savefig('fig/vx_out_integral.png')

	aspect=1.0
	perc=93
	arr = vz_out_integral.T
	nt,nr = arr.shape
	ratio = float(nr)/nt
	imshow_aspect = ratio/aspect
	plt.figure()
	im = plt.imshow(arr, interpolation='nearest', cmap=cmap, aspect=imshow_aspect, vmax=vmax, vmin=vmin)  
	plt.colorbar()
	plt.title('vz_conv')
	plt.xlabel('nx')
	plt.ylabel('nt')
	plt.savefig('fig/vz_out_integral.png')
 
	#########################################
	######################################### 
	
def plot_traces(shotgather_extra, geom, shotgather_wf, xp, time_i, time_f):
	###load extrapolate wavefield
	dict_extrapolate =spio.loadmat(shotgather_extra)
	vx_ex   = dict_extrapolate['vx_extrapolate']
	vz_ex   = dict_extrapolate['vz_extrapolate']
	vx_out_integral  = dict_extrapolate['vx_out_integral']
	vz_out_integral  = dict_extrapolate['vz_out_integral']
	##load source
	dict_geom = spio.loadmat(geom)
	source_data  = dict_geom['source_data'][0]
	dt           = dict_geom['dt'][0][0]
	## load shotgather_wf
	dict_geom = spio.loadmat(shotgather_wf)
	vx_wf     = dict_geom['vx_wf']
	vz_wf     = dict_geom['vz_wf']
	nt, nr = vx_wf.shape
	##define traces position to plot
	trace_x_position = xp
	
	##VX FULL
	plt.figure()
	im=plt.plot(vx_wf[:,xp],'g', label='Full',     linewidth=2) 
	plt.plot   (vx_ex[:,xp],'r', label='Local',  linewidth=2) 
	plt.title('vx, trace=%.2f'%(xp))
	plt.ylabel('Amplitude' )
	plt.legend()
	plt.xlabel('nt')
	plt.savefig('fig/traces/vx_trace=%.2f.png'%(xp))
	
	##VZ FULL 
	plt.figure()
	im=plt.plot(vz_wf[:,xp],'g', label='Full',     linewidth=2) 
	plt.plot   (vz_ex[:,xp],'r', label='Local',  linewidth=2) 
	plt.title('vz, trace=%.2f'%(xp))
	plt.ylabel('Amplitude' )
	plt.legend()
	plt.xlabel('nt')
	plt.savefig('fig/traces/vz_trace=%.2f.png'%(xp))
	
	#plot reflection of interest
	
	itime=time_i# time initial
	itimef=time_f# time final
	##VX
	xx = np.arange(itime,itimef,1)*dt/1000
	yx1 = vx_wf[itime:itimef,xp]
	yx2 = vx_ex[itime:itimef,xp]
	
	plt.figure()
	im=plt.plot(xx, yx1,'k',   label='Full',   linewidth=1)# 
	plt.plot   (xx, yx2,'*--r',  label='Local',  linewidth=1)# 
	plt.ylabel('Amplitude')
	plt.xlabel('Time (s)')
	plt.legend()
	plt.tight_layout()
	plt.savefig('fig/traces/vx_short_trace=%.2f.png'%(xp))
	
	##VZ
	xz  = np.arange(itime,itimef,1)*dt/1000
	yz1 = vz_wf[itime:itimef,xp]
	yz2 = vz_ex[itime:itimef,xp]
	
	plt.figure()
	im=plt.plot(xz, yz1,'k',   label='Full',   linewidth=1)# 
	plt.plot   (xz, yz2,'*--r',  label='Local',  linewidth=1)# 
	plt.ylabel('Amplitude' )
	plt.legend()
	plt.xlabel('Time (s)')
	plt.tight_layout()
	plt.savefig('fig/traces/vz_short_trace=%.2f.png'%(xp))
	
	
def simple_model_plot(arr1, arr2, arr3, name_output):
	cmap = plt.get_cmap('seismic')
	plt.figure()
	plt.subplot(311)
	im = plt.imshow(arr1, cmap=cmap, aspect='auto')
	plt.colorbar()
	
	plt.subplot(312)
	im = plt.imshow(arr2, cmap=cmap, aspect='auto')
	plt.colorbar()
	
	plt.subplot(313)
	im = plt.imshow(arr3, cmap=cmap, aspect='auto')
	plt.colorbar()	
		
	plt.savefig('fig/%s'%name_output)
	
	
def plot_simple_shotgaher(arr, name_output):
	plt.figure()
	im=plt.imshow(arr, cmap='seismic', aspect='auto', vmax=0.02, vmin=-0.02)
	plt.colorbar()	
	plt.savefig('fig/%s'%name_output)
	
def plot_simple_traces(arr_full, arr_local, xp, time_i, time_f ):
	
	##load full
	dict_full = spio.loadmat(arr_full)
	vx_full     = dict_full['vx']

	###load local
	dict_local =spio.loadmat(arr_local)
	vx_local  = dict_local['vx']
	
	##VX FULL
	plt.figure()
	im=plt.plot(vx_full [time_i:time_f,xp], 'k', label='Full',   linewidth=2) 
	plt.plot   (vx_local[time_i+1:time_f+1,xp], '*-r', label='Local',  linewidth=2) 
	plt.ylabel('Amplitude' )
	plt.legend()
	plt.xlabel('time steps')
	plt.savefig('fig/trace=%.2f.png'%(xp))
	
	
	
	



