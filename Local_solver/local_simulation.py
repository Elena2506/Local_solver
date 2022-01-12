"""
    Local Solver

    A Python program to inject wavefield as point sources
    
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
__all__ = ['local_simulation']
tracemalloc.start()


###############################################################################
## Load Point sources to inject
###############################################################################


def local_simulation(src_coordinate, so, model_dirt, mi):

    """
    Local simulation 

    Parameters:

    * src_coordinate: 1D array source coordinate
        
    * so: space-order   
    
    * model_dirt: model directory address
    
    * mi: non perturbed (nonp) or perturbed (pert) model 
    
 
    Returns:

    * Extrapolation surface wavefields
    * v1 data for plotting

    """
    ################################################################################
    ##Load local model
    #################################################################################
    dict_model =spio.loadmat( model_dirt)
    vp      = dict_model['vp']
    vs      = dict_model['vs']
    rho     = dict_model['rho']
    space   = dict_model['spacing']
    nbpml   = dict_model['nbl'][0][0]
    t0      = dict_model['t0'][0][0]
    tn      = dict_model['tn'][0][0]
    origin  = (0., 0.) 
    nbl     = nbpml
    spacing = (space[0][0], space[0][1])
    shape   = (vp.shape[0], vp.shape[1])
    dx = space[0][0]
    dz = space[0][1]
    
    model = ModelElastic(vp=vp,vs=vs, b=1./rho, origin=origin, shape=shape, spacing=spacing, space_order=so, nbl=nbl)
    
    ###################################
    #Time
    ###################################
    dt          = model.critical_dt
    time_range  = TimeAxis(start=t0, stop=tn, step=dt)
    save        = time_range.num

    ###############################################################################
    ## Load injection surface (S) coordinates
    ###############################################################################
    
    dict_inj =spio.loadmat('output/coord_inj_local.mat') 
    inj_coordinates_t = dict_inj['inj_coordinates_t']
    inj_coordinates_l = dict_inj['inj_coordinates_l']
    inj_coordinates_b = dict_inj['inj_coordinates_b']
    inj_coordinates_r = dict_inj['inj_coordinates_r']
    inj_t            =  dict_inj['inj_t']
    nx_inj           =  dict_inj['nx_inj'][0][0]
    nz_inj           =  dict_inj['nz_inj'][0][0]
    
    
    for i in src_coordinate:
    	print ('injection for a source located at x=', i[0])
    	dict_pointsource =spio.loadmat('output/point_sources_src_nonp_'+str(i[0])+'0.mat')
    	vz_inj_t  = dict_pointsource['vz_inj_t']
    	vx_inj_t  = dict_pointsource['vx_inj_t']
    	txx_inj_t = dict_pointsource['txx_inj_t']
    	txz_inj_t = dict_pointsource['txz_inj_t']
    	tzz_inj_t = dict_pointsource['tzz_inj_t']
    	
    	vz_inj_l  = dict_pointsource['vz_inj_l']
    	vx_inj_l  = dict_pointsource['vx_inj_l']
    	txx_inj_l = dict_pointsource['txx_inj_l']
    	txz_inj_l = dict_pointsource['txz_inj_l']
    	tzz_inj_l = dict_pointsource['tzz_inj_l']
    	
    	vz_inj_b  = dict_pointsource['vz_inj_b']
    	vx_inj_b  = dict_pointsource['vx_inj_b']
    	txx_inj_b = dict_pointsource['txx_inj_b']
    	txz_inj_b = dict_pointsource['txz_inj_b']
    	tzz_inj_b = dict_pointsource['tzz_inj_b']
    	
    	vz_inj_r  = dict_pointsource['vz_inj_r']
    	vx_inj_r  = dict_pointsource['vx_inj_r']
    	txx_inj_r = dict_pointsource['txx_inj_r']
    	txz_inj_r = dict_pointsource['txz_inj_r']
    	tzz_inj_r = dict_pointsource['tzz_inj_r']
    	
    	
    	###################################
    	#Forward modeling
    	###################################
    	# Now we create the velocity and pressure fields
    	
    	x, z = model.grid.dimensions
    	t = model.grid.stepping_dim
    	time = model.grid.time_dim
    	s = time.spacing
    	
    	v   = VectorTimeFunction(name='v', grid=model.grid, space_order=so, time_order=2, save=time_range.num)
    	tau = TensorTimeFunction(name='t', grid=model.grid, space_order=so, time_order=2, save=time_range.num)
    	
    	# Now let's try and create the staggered updates
    	# Lame parameters
    	l, mu, ro = model.lam, model.mu, model.b
    	#################################################################
    	
    	# # # Equations and Stencils
    	u_v = Eq(v.forward,    model.damp * (v + s*ro*div(tau)))
    	u_t = Eq(tau.forward,  model.damp *  (tau + s * (l * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).T))))
    	
    	##injection top
    	################################################################################
    	source_vx    = PointSource(name='source_vx', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_t, data=vx_inj_t)
    	src_in_vx    = source_vx.inject(field=v.forward[0], expr=source_vx)
    	
    	source_vz    = PointSource(name='source_vz', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_t, data=vz_inj_t)
    	src_in_vz    = source_vz.inject(field=v.forward[1], expr=source_vz)
    	
    	source_txx   = PointSource(name='source_txx', grid=model.grid, time_range=time_range,coordinates=inj_coordinates_t, data=txx_inj_t)
    	src_in_txx   = source_txx.inject(field=tau.forward[0, 0], expr=source_txx)
    	
    	source_tzz   = PointSource(name='source_tzz', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_t, data=tzz_inj_t)
    	src_in_tzz   = source_tzz.inject(field=tau.forward[1, 1], expr=source_tzz)
    	
    	source_txz   = PointSource(name='source_txz', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_t, data=txz_inj_t)
    	src_in_txz   = source_txz.inject(field=tau.forward[0, 1], expr=source_txz)
    	
    	##injection left
    	
    	################################################################################
    	source_vx_l    = PointSource(name='source_vx_l', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_l, data=vx_inj_l)
    	src_in_vx_l    = source_vx_l.inject(field=v.forward[0], expr=source_vx_l)
    	
    	source_vz_l    = PointSource(name='source_vz_l', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_l, data=vz_inj_l)
    	src_in_vz_l    = source_vz_l.inject(field=v.forward[1], expr=source_vz_l)
    	
    	source_txx_l   = PointSource(name='source_txx_l', grid=model.grid, time_range=time_range,coordinates=inj_coordinates_l, data=txx_inj_l)
    	src_in_txx_l   = source_txx_l.inject(field=tau.forward[0, 0], expr=source_txx_l)
    	
    	source_tzz_l   = PointSource(name='source_tzz_l', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_l, data=tzz_inj_l)
    	src_in_tzz_l   = source_tzz_l.inject(field=tau.forward[1, 1], expr=source_tzz_l)
    	
    	source_txz_l   = PointSource(name='source_txz_l', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_l, data=txz_inj_l)
    	src_in_txz_l   = source_txz_l.inject(field=tau.forward[0, 1], expr=source_txz_l)
    	
    	##injection bottom
    	################################################################################
    	
    	source_vx_b    = PointSource(name='source_vx_b', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_b, data=vx_inj_b)
    	src_in_vx_b    = source_vx_b.inject(field=v.forward[0], expr=source_vx_b)
    	
    	source_vz_b    = PointSource(name='source_vz_b', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_b, data=vz_inj_b)
    	src_in_vz_b    = source_vz_b.inject(field=v.forward[1], expr=source_vz_b)
    	
    	source_txx_b   = PointSource(name='source_txx_b', grid=model.grid, time_range=time_range,coordinates=inj_coordinates_b, data=txx_inj_b)
    	src_in_txx_b   = source_txx_b.inject(field=tau.forward[0, 0], expr=source_txx_b)
    	
    	source_tzz_b   = PointSource(name='source_tzz_b', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_b, data=tzz_inj_b)
    	src_in_tzz_b   = source_tzz_b.inject(field=tau.forward[1, 1], expr=source_tzz_b)
    	
    	source_txz_b   = PointSource(name='source_txz_b', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_b, data=txz_inj_b)
    	src_in_txz_b   = source_txz_b.inject(field=tau.forward[0, 1], expr=source_txz_b)
    	
    	##injection right
    	################################################################################
    	
    	source_vx_r    = PointSource(name='source_vx_r', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_r, data=vx_inj_r)
    	src_in_vx_r    = source_vx_r.inject(field=v.forward[0], expr=source_vx_r)
    	
    	source_vz_r    = PointSource(name='source_vz_r', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_r, data=vz_inj_r)
    	src_in_vz_r    = source_vz_r.inject(field=v.forward[1], expr=source_vz_r)
    	
    	source_txx_r   = PointSource(name='source_txx_r', grid=model.grid, time_range=time_range,coordinates=inj_coordinates_r, data=txx_inj_r)
    	src_in_txx_r   = source_txx_r.inject(field=tau.forward[0, 0], expr=source_txx_r)
    	
    	source_tzz_r   = PointSource(name='source_tzz_r', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_r, data=tzz_inj_r)
    	src_in_tzz_r   = source_tzz_r.inject(field=tau.forward[1, 1], expr=source_tzz_r)
    	
    	source_txz_r   = PointSource(name='source_txz_r', grid=model.grid, time_range=time_range, coordinates=inj_coordinates_r, data=txz_inj_r)
    	src_in_txz_r   = source_txz_r.inject(field=tau.forward[0, 1], expr=source_txz_r)
    	
    	src_in_expre = (src_in_vx   + src_in_vz   + src_in_txx   + src_in_tzz   + src_in_txz+
    					src_in_vx_l + src_in_vz_l + src_in_txx_l + src_in_tzz_l + src_in_txz_l +
    					src_in_vx_b + src_in_vz_b + src_in_txx_b + src_in_tzz_b + src_in_txz_b + 
    					src_in_vx_r + src_in_vz_r + src_in_txx_r + src_in_tzz_r + src_in_txz_r)
    					
    	op = Operator([u_v] + [u_t] + src_in_expre, save=True)
    	op(dt=model.critical_dt)
    	
    	###############################################################################
    	### Load extrapolation surface location (Se)
    	###############################################################################
    	
    	dict_se=spio.loadmat('output/coord_se_local.mat') 
    	se_coordinates = dict_se['se_coord']
    	n_se           = dict_se['n_se'][0][0]
    	
    	#####################################
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
    		
    	##Save wavefield at Se local
    	outdict = dict()
    	outdict['tzz']= tzz_se
    	outdict['txx']= txx_se
    	outdict['txz']= txz_se
    	outdict['vz']= vz_se
    	outdict['vx']= vx_se
    	outdict['tn']= tn
    	spio.savemat('output/waveform_se_local_src_%s_%.2f.mat'%(mi,i[0]), outdict)
    	
    	outdict = dict()
    	outdict['v1']= v[1].data[int(save/2),:,:]
    	spio.savemat('output/v1_local_%s.mat'%mi, outdict)
    	#memory 
    	current, peak = tracemalloc.get_traced_memory()
    	print(f"Current memory usage in the injection forward modeling is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    	tracemalloc.stop()
    	
    	print('Finish local simulation')
