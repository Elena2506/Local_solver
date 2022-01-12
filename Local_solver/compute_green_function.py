"""
    Local Solver

    A Python program to compute Green functions 
    
    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""

import scipy.io as spio
from devito import *
from examples.seismic import *

__all__ = ['Green_function']

def Green_function(model, geom, so, direction, save_dir):

    """
    Use the forward modeling to compute the Green functions between receiver positions and Se

    Parameters:

    * model: address to input model
        
    * geom: address to input geometry information 

    * so: space-order   
    
    * direction: x or z direction to inject 
    
    * save_dir: location of the save directory 
    
 
    Returns:

    * Green function library

    """
    
    ##################################
    ###LOAD MODEL
    ##################################
    dict_fullmodel =spio.loadmat(model)
    
    vp      = dict_fullmodel['vp']
    vs      = dict_fullmodel['vs']
    rho     = dict_fullmodel['rho']
    space   = dict_fullmodel['spacing']
    shaping = dict_fullmodel['shape']  
    nbl     = dict_fullmodel['nbl'][0][0]
    t0      = dict_fullmodel['t0'][0][0]
    tn      = dict_fullmodel['tn'][0][0] 
    origin  = (0., 0.) 
    nbl=100
    nbpml   = nbl
    
    spacing = (space[0][0], space[0][1])
    shape   = (shaping[0][0], shaping[0][1])
    dx=space[0][0]
    dz=space[0][1]
    
    ################################################
    dict_geom = spio.loadmat(geom)
    ################################################
    rec_coord      =dict_geom['rec_coord']
    se_coordinates =dict_geom['se_coord']
    f0             =dict_geom['f0']
    ###############################################
    
    ##################################
    ###Define Model for Devito 
    ##################################
    so=so
    model = ModelElastic(vp=vp,vs=vs,b=1./rho, origin=origin, shape=shape, spacing=spacing, space_order=so, nbl=nbl)
    
    ##################################
    ###Define geometry
    ##################################
    
    ################### time
    dt = model.critical_dt
    time_range = TimeAxis(start=t0, stop=tn, step=dt)
    save=time_range.num
    ############################################################
    
    #############Forward modeling 
    for i in rec_coord:
    	print('Source',i)
    	src = RickerSource(name='src', grid=model.grid, f0=f0, time_range=time_range)
    	src.coordinates.data[:] = i
    	
    	x_src = int(i[0]/dx)
    	z_src = int(i[1]/dz)
    	rho_src_position  = dict_fullmodel['rho'][x_src][z_src]
    	
    	src.data[:] = src.data/(rho_src_position)
    	
    	x, z = model.grid.dimensions
    	t = model.grid.stepping_dim
    	time = model.grid.time_dim
    	s = time.spacing
    	v   = VectorTimeFunction(name='v', grid=model.grid, space_order=so, time_order=2, save=time_range.num)
    	tau = TensorTimeFunction(name='t', grid=model.grid, space_order=so, time_order=2, save=time_range.num)
    	
    	src_dir = src.inject(field=v.forward[direction], expr=s*src)
    	
    	
    	#Lame parameters
    	l, mu, ro = model.lam, model.mu, model.b
    	
    	#fdelmodc reference implementation
    	u_v = Eq(v.forward,    model.damp * (v + s*ro*div(tau)))
    	u_t = Eq(tau.forward,  model.damp *  (tau + s * (l * diag(div(v.forward)) + mu * (grad(v.forward) + grad(v.forward).T))))
    	
    	op = Operator([u_v] + [u_t] + src_dir ,  save=True)
    	op(dt=model.critical_dt)
    	
    	##################################################################################
    	####### extrapolation surface position
    	####### record the wavefield at the extrapolation surface for each source at  
    	####### receiver position. The source is injected in x or z direction. 
    	##################################################################################
    	
    	txx_se = np.zeros((save, se_coordinates.shape[0]))
    	txz_se = np.zeros((save, se_coordinates.shape[0]))
    	tzz_se = np.zeros((save, se_coordinates.shape[0]))
    	vx_se  = np.zeros((save, se_coordinates.shape[0]))
    	vz_se  = np.zeros((save, se_coordinates.shape[0]))
    	
    	se_len=len(se_coordinates)
    	
    	for r in range (se_len):
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
    	
    	isorc= int(src.coordinates.data[0][0])
    	outdict = dict()
    	outdict['tzz']= tzz_se
    	outdict['txx']= txx_se
    	outdict['txz']= txz_se
    	outdict['vz']= vz_se
    	outdict['vx']= vx_se 
    	spio.savemat('output/G_V%s/Component_source_greens_se%.2f'%(save_dir,isorc), outdict)
    	
    	print('Finished Green function direction %s'%save_dir)
		##################################################################################




