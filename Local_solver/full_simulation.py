"""
    Local Solver

    A Python program to create forward modeling in the full domain
    
    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""

import numpy as np
import scipy.io as spio
from devito import *
from examples.seismic import *
import tracemalloc

__all__ = ['full_simulation']

tracemalloc.start()

"""
	Full simulation 

    Parameters:

    * src_coordinate: 1D array source coordinate
        
    * so: space-order   
    
    * model_dirt: model directory address
    
    * mi: non perturbed (nonp) or perturbed (pert) model 
    
 
    Returns:
    
    * geometry
    * Point Sources
    * Shotgather
    * Extrapolation surface wavefields
    
    """
def full_simulation(src_coordinate, so, model_dirt, mi):
    dict_fullmodel =spio.loadmat(model_dirt)
    
    vp      = dict_fullmodel['vp']
    vs      = dict_fullmodel['vs']
    rho     = dict_fullmodel['rho']
    space   = dict_fullmodel['spacing']
    nbl     = dict_fullmodel['nbl'][0][0]
    t0      = dict_fullmodel['t0'][0][0]
    tn      = dict_fullmodel['tn'][0][0]
    origin  = (0., 0.) 
    nbpml=nbl
    shape   = (vp.shape[0],vp.shape[1])
    spacing = (space[0][0], space[0][1])
    dx = space[0][0]
    dz = space[0][1]
    
    model = ModelElastic(vp=vp, vs=vs, b=1./rho, origin=origin, shape=shape, spacing=spacing, space_order=so, nbl=nbl)
    
    ###################################
    #Time
    ###################################
    
    dt = model.critical_dt
    time_range = TimeAxis(start=t0, stop=tn, step=dt)
    save=time_range.num
    
    for i in src_coordinate:
    	f0=0.015
    	src = RickerSource(name='src', grid=model.grid, f0=f0, time_range=time_range)
    	src.coordinates.data[:]  = i
    	print('source at', src.coordinates.data)
    	source_data = np.asarray(src.data)
    	source_data = np.reshape(source_data,source_data.size)
    	#####################################################################################
    	
    	###################################
    	#Forward modeling
    	###################################
    	# Now we create the velocity and pressure fields
    	clear_cache()
    	
    	x, z = model.grid.dimensions
    	t = model.grid.stepping_dim
    	time = model.grid.time_dim
    	s = time.spacing
    	v   = VectorTimeFunction(name='v', grid=model.grid, space_order=so, time_order=2, save=time_range.num)
    	tau = TensorTimeFunction(name='t', grid=model.grid, space_order=so, time_order=2, save=time_range.num)
    	
    	# The source injection term
    	src_xx = src.inject(field=tau.forward[0, 0], expr=s*src)
    	src_zz = src.inject(field=tau.forward[1, 1], expr=s*src)
    	# Now let's try and create the staggered updates
    	# Lame parameters
    	l, mu, ro = model.lam, model.mu, model.b
    	# fdelmodc reference implementation
    	u_v = Eq(v.forward,    model.damp * (v + s*ro*div(tau)))
    	u_t = Eq(tau.forward,  model.damp *  (tau + s * (l * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).T))))
    	
    	op = Operator([u_v] + [u_t] + src_xx + src_zz, save=True)
    	
    	## using the dt critical and call for the Operator devito/operator 
    	
    	op(dt=model.critical_dt)
    	
    	##################################################
    	###Shotgather
    	##################################################
    	
    	##load receiver locations
    	dict_rec   = spio.loadmat('output/coord_rec.mat') 
    	n_rec     = dict_rec['nrec'][0][0]
    	xmin_rec  = dict_rec['xmin_rec'][0][0]
    	xmax_rec  = dict_rec['xmax_rec'][0][0]
    	drec      = dict_rec['drec'][0][0]
    	z_rec     = dict_rec['z_rec'][0][0]
    	
    	rec = Receiver(name="rec", grid=model.grid, npoint=int(n_rec), time_range=time_range)
    	rec.coordinates.data[:, 0]  = np.arange(xmin_rec, xmax_rec, drec)
    	rec.coordinates.data[:, -1] = z_rec
    	
    	#################################################
    	vx_rec  = np.zeros((save, n_rec))
    	vz_rec  = np.zeros((save, n_rec))
    	
    	for i in range (n_rec):
    		xrec=int((rec.coordinates_data[i][0])/dx+nbl)
    		zrec=int((rec.coordinates_data[i][1])/dz+nbl)
    		vx_t2 = v[0].data[:,xrec,zrec]
    		vx_rec[:,i]=vx_t2
    		vz_t2 = v[1].data[:,xrec,zrec]
    		vz_rec[:,i]=vz_t2
    	##Save shotgather
    	isrc = int(src.coordinates.data[0][0])
    	outdict       = dict()
    	outdict['vz'] = vz_rec
    	outdict['vx'] = vx_rec
    	spio.savemat('output/shotgather_src_%s_%.2f.mat'%(mi, isrc),  outdict)
    	
    	##################################################
    	### Extrapolation surface location (Se)
    	##################################################
    	
    	##load coordinates extrapolation surface Se
    	dict_se =spio.loadmat('output/coord_se.mat')
    	se_coordinates = dict_se['se_coord']
    	n_se           =dict_se['n_se'][0][0]
    	
    	txx_se = np.zeros((save, n_se))
    	txz_se = np.zeros((save, n_se))
    	tzz_se = np.zeros((save, n_se))
    	vx_se  = np.zeros((save, n_se))
    	vz_se  = np.zeros((save, n_se))
    	
    	for r in range (n_se):
    		xse=int((se_coordinates[r][0])/dx+nbpml)
    		zse=int((se_coordinates[r][1])/dz+nbpml)
    		txx_se_i = tau[0,0].data[:,xse,zse]
    		txx_se[:,r] = txx_se_i
    		tzz_se_i = tau[1,1].data[:,xse,zse]
    		tzz_se[:,r] = tzz_se_i
    		vx_se_i = v[0].data[:,xse,zse]
    		vx_se[:,r]=vx_se_i
    		txz_se_i = tau[0,1].data[:,xse,zse]
    		txz_se[:,r] = txz_se_i
    		vz_se_i = v[1].data[:,xse,zse]
    		vz_se[:,r]=vz_se_i
    	
    	##Save wavefield at Se
    	outdict = dict()
    	outdict['tzz']= tzz_se
    	outdict['txx']= txx_se
    	outdict['txz']= txz_se
    	outdict['vz']= vz_se
    	outdict['vx']= vx_se
    	outdict['tn']= tn
    	spio.savemat('output/waveform_se_full_src_%s_%.2f.mat'%(mi, isrc), outdict)
    	
    	###############################################################################
    	## Injection surface (S) 
    	## For each part of the injection surface (i.e. top, left, bottom and right) compute:
    	## 1. Take wavefields in the injection surface (S) positions 
    	## 2. extract the values of lambda, 1./rho and mu 
    	## 3. Use equation (2) to compute the point-sources
    	###############################################################################
    	
    	##load coordinates injection surface (S)
    	dict_inj =spio.loadmat('output/coord_inj.mat') 
    	inj_coordinates_t =dict_inj['inj_coordinates_t']
    	inj_coordinates_l =dict_inj['inj_coordinates_l']
    	inj_coordinates_b =dict_inj['inj_coordinates_b']
    	inj_coordinates_r =dict_inj['inj_coordinates_r']
    	inj_t  = dict_inj['inj_t']
    	nx_inj = dict_inj['nx_inj'][0][0]
    	nz_inj = dict_inj['nz_inj'][0][0]
    	
    	##Compute Source - injection## 
    	print('Computing new point sources to inject') 
    	
    	###############################################################################   
    	## Top
    	###############################################################################
    	
    	txx_inj_t = np.zeros((save, nx_inj))
    	txz_inj_t = np.zeros((save, nx_inj))
    	tzz_inj_t = np.zeros((save, nx_inj))
    	vx_inj_t  = np.zeros((save, nx_inj))
    	vz_inj_t  = np.zeros((save, nx_inj))
    	
    	for r in range (nx_inj):
    		xinj=int((inj_coordinates_t[r][0])/dx+nbpml)
    		zinj=int((inj_coordinates_t[r][1])/dz+nbpml)
    		
    		txx_t = tau[0,0].data[:,xinj,zinj]
    		tzz_t = tau[1,1].data[:,xinj,zinj]
    		txz_t = tau[0,1].data[:,xinj,zinj]
    		vx_t = v[0].data[:,xinj,zinj]
    		vz_t = v[1].data[:,xinj,zinj]
    		b_top    =ro.data[xinj,zinj]
    		la_top   =l.data [xinj,zinj]
    		mum_top  =mu.data[xinj,zinj]
    		vx_inj_t[:,r]=b_top*(dt/dz)*(txz_t)*(-1)
    		vz_inj_t[:,r]=b_top*(dt/dz)*(tzz_t)*(-1)
    		txx_inj_t[:,r]=(dt/dz)*(la_top*vz_t)*(-1) 
    		txz_inj_t[:,r]=(dt/dz)*(mum_top*vx_t)*(-1)
    		tzz_inj_t[:,r]=(dt/dz)*(la_top+2*mum_top)*vz_t*(-1)
    		
    	###############################################################################
    	##Left
    	###############################################################################
    	
    	txx_inj_l = np.zeros((save, nz_inj))
    	txz_inj_l = np.zeros((save, nz_inj))
    	tzz_inj_l = np.zeros((save, nz_inj))
    	vx_inj_l  = np.zeros((save, nz_inj))
    	vz_inj_l  = np.zeros((save, nz_inj))
    	
    	for r in range (nz_inj):
    		xinj=int((inj_coordinates_l[r][0])/dx+nbpml)
    		zinj=int((inj_coordinates_l[r][1])/dz+nbpml)
    		txx_le = tau[0,0].data[:,xinj,zinj]
    		tzz_le = tau[1,1].data[:,xinj,zinj]
    		txz_le = tau[0,1].data[:,xinj,zinj]
    		vx_le = v[0].data[:,xinj,zinj]
    		vz_le = v[1].data[:,xinj,zinj]
    		b_left   =ro.data[xinj,zinj]
    		la_left  =l.data [xinj,zinj]
    		mum_left =mu.data[xinj,zinj]
    		
    		vx_inj_l[:,r]=b_left*(dt/dz)*(txx_le)*(-1)
    		vz_inj_l[:,r]=b_left*(dt/dz)*(txz_le)*(-1)
    		tzz_inj_l[:,r]=(dt/dz)*la_left*vx_le*(-1)
    		txz_inj_l[:,r]=(dt/dz)*mum_left*vz_le*(-1)
    		txx_inj_l[:,r]=(dt/dz)*(la_left+2*mum_left)*vx_le*(-1)
    		
    	###############################################################################
    	##Bottom 
    	###############################################################################
    	
    	txx_inj_b = np.zeros((save, nx_inj))
    	txz_inj_b = np.zeros((save, nx_inj))
    	tzz_inj_b = np.zeros((save, nx_inj))
    	vx_inj_b  = np.zeros((save, nx_inj))
    	vz_inj_b  = np.zeros((save, nx_inj))
    	
    	for r in range (nx_inj):
    		xinj=int((inj_coordinates_b[r][0])/dx+nbpml)
    		zinj=int((inj_coordinates_b[r][1])/dz+nbpml)
    		
    		txx_b = tau[0,0].data[:,xinj,zinj]
    		tzz_b = tau[1,1].data[:,xinj,zinj]
    		txz_b = tau[0,1].data[:,xinj,zinj]
    		vx_b = v[0].data[:,xinj,zinj]
    		vz_b = v[1].data[:,xinj,zinj]
    		b_b=ro.data[xinj,zinj]
    		la_b=l.data[xinj,zinj]
    		mum_b=mu.data[xinj,zinj]
    		vx_inj_b[:,r]=b_b*(dt/dz)*(txz_b)
    		vz_inj_b[:,r]=b_b*(dt/dz)*(tzz_b)
    		txx_inj_b[:,r]=(dt/dz)*la_b*vz_b
    		txz_inj_b[:,r]=(dt/dz)*mum_b*vx_b
    		tzz_inj_b[:,r]=(dt/dz)*(la_b+2*mum_b)*vz_b
    		
    	###############################################################################
    	##Rigth 
    	###############################################################################
    	
    	txx_inj_r = np.zeros((save, nz_inj))
    	txz_inj_r = np.zeros((save, nz_inj))
    	tzz_inj_r = np.zeros((save, nz_inj))
    	vx_inj_r  = np.zeros((save, nz_inj))
    	vz_inj_r  = np.zeros((save, nz_inj))
    	
    	for r in range (nz_inj):
    		xinj=int((inj_coordinates_r[r][0])/dx+nbpml)
    		zinj=int((inj_coordinates_r[r][1])/dz+nbpml)
    		
    		txx_re = tau[0,0].data[:,xinj,zinj]
    		tzz_re = tau[1,1].data[:,xinj,zinj]
    		txz_re = tau[0,1].data[:,xinj,zinj]
    		vx_re = v[0].data[:,xinj,zinj]
    		vz_re = v[1].data[:,xinj,zinj]
    		b_r   = ro.data[xinj,zinj]
    		la_r  =  l.data[xinj,zinj]
    		mum_r = mu.data[xinj,zinj]
    		
    		vx_inj_r[:,r]  = b_r*(dt/dz)*(txx_re)
    		vz_inj_r[:,r]  = b_r*(dt/dz)*(txz_re)
    		tzz_inj_r[:,r] = (dt/dz)*la_r*vx_re
    		txz_inj_r[:,r]  =(dt/dz)*mum_r*vz_re
    		txx_inj_r[:,r] = (dt/dz)*(la_r+2*mum_r)*vx_re
    		
    	print('finish the point sources calculations') 
    	###############################################################
    	##save point sources
    	
    	outdict = dict()
    	outdict['vz_inj_t']   = vz_inj_t
    	outdict['vx_inj_t']   = vx_inj_t
    	outdict['txx_inj_t']  = txx_inj_t
    	outdict['txz_inj_t']  = txz_inj_t
    	outdict['tzz_inj_t']  = tzz_inj_t
    	outdict['vz_inj_l']   = vz_inj_l
    	outdict['vx_inj_l']   = vx_inj_l
    	outdict['txx_inj_l']  = txx_inj_l
    	outdict['txz_inj_l']  = txz_inj_l
    	outdict['tzz_inj_l']  = tzz_inj_l
    	outdict['vz_inj_b']   = vz_inj_b
    	outdict['vx_inj_b']   = vx_inj_b
    	outdict['txx_inj_b']  = txx_inj_b
    	outdict['txz_inj_b']  = txz_inj_b
    	outdict['tzz_inj_b']  = tzz_inj_b
    	outdict['vz_inj_r']   = vz_inj_r
    	outdict['vx_inj_r']   = vx_inj_r
    	outdict['txx_inj_r']  = txx_inj_r
    	outdict['txz_inj_r']  = txz_inj_r
    	outdict['tzz_inj_r']  = tzz_inj_r
    	spio.savemat('output/point_sources_src_%s_%.2f.mat'%(mi, isrc), outdict)
    	
    	print('saved point injection sources') 
    	###################################################################################
    	## save geometry parameters 
    	outdict                 = dict()
    	outdict['nbl']          = nbl
    	outdict['src_coord']    = src.coordinates_data
    	outdict['rec_coord']    = rec.coordinates.data
    	outdict['se_coord']     = se_coordinates
    	outdict['spacing']      = spacing
    	outdict['dt']           = dt
    	outdict['f0']           = f0
    	outdict['source_data']  = source_data
    	outdict['v1']           = v[1].data[int(save/2),:,:]
    	spio.savemat('output/geometry_%s_%.2f.mat'%(mi, isrc), outdict)
    	
    	#memory for each shot
    	current, peak = tracemalloc.get_traced_memory()
    	print(f"Current memory usage during forward modeling is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    	tracemalloc.stop()  
    	print('Finished full simulation')
		
