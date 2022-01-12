"""
    Local Solver

    A Python program to:
    - Define elastic velocity model from SEAM 2D data.
    - Define Full domain.  

    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""


import numpy as np
import scipy.io as spio
from devito import *
from examples.seismic import *

__all__ = ['full_model', 'full_model_pert', 'local_model']



def full_model(nbl, t0, tn, space, SM_vp, SM_vs, SM_rho, xmin, xmax, zmin, zmax, xmin_inj, xmax_inj, zmin_inj, zmax_inj, xmin_se, xmax_se, zmin_se,
				reflector_depth, xmax_local, xmin_local, zmin_local, zmax_local, xmin_inj_local, xmax_inj_local, zmin_inj_local, zmax_inj_local, 
				xmax_se_local,  xmin_se_local, zmin_se_local):
	
	vp_max =np.max(SM_vp)
	vs_max =np.max(SM_vs)
	rho_max =np.max(SM_rho)
	vp = SM_vp
	vp[0:1, 0:1]= vp_max
	vs = SM_vs
	vs[0:1, 0:1]= vs_max
	rho = SM_rho
	rho[0:1, 0:1]= rho_max
	
	shape   = (vp.shape[0],vp.shape[1])
	# print(np.max(vp), np.max(vs), np.max(rho))
	##################################################################################
	#Save full model 
	###################################################################################
	 
	outdict = dict()
	outdict['vp']     = vp
	outdict['vs']     = vs
	outdict['rho']    = rho
	outdict['spacing']= space
	outdict['shape']  = vp.shape
	outdict['nbl']    = nbl
	outdict['t0']     = t0
	outdict['tn']     = tn
	##Local domain index in the full domain
	outdict['xmin']=xmin
	outdict['xmax']=xmax
	outdict['zmin']=zmin
	outdict['zmax']=zmax
	##Injection full index
	outdict['xmin_inj']=xmin_inj
	outdict['xmax_inj']=xmax_inj
	outdict['zmin_inj']=zmin_inj
	outdict['zmax_inj']=zmax_inj
	##Extrapolation surface index
	outdict['xmin_se']   =xmin_se
	outdict['xmax_se']   =xmax_se
	outdict['zmin_se']   =zmin_se
	## Reflector position
	outdict['reflector_depth']   = reflector_depth
	## Local domain index
	outdict['xmax_local']=xmax_local
	outdict['xmin_local']=xmin_local
	outdict['zmin_local']=zmin_local
	outdict['zmax_local']=zmax_local
	##injection local index
	outdict['xmin_inj_local']=xmin_inj_local
	outdict['xmax_inj_local']=xmax_inj_local
	outdict['zmin_inj_local']=zmin_inj_local
	outdict['zmax_inj_local']=zmax_inj_local
	##extrapolation surface local index
	outdict['xmax_se_local']=xmax_se_local
	outdict['xmin_se_local']=xmin_se_local
	outdict['zmin_se_local']=zmin_se_local
	# Reflector position
	outdict['reflector_depth_local']   = reflector_depth_local
	
	spio.savemat('output/full.mat', outdict)
	return vp, vs, rho, shape, spacing
	
#######################

def full_model_pert(nbl, t0, tn, space, SM_vp, SM_vs, SM_rho, xmin, xmax, zmin, zmax, xmin_inj, xmax_inj, zmin_inj, zmax_inj, xmin_se, xmax_se, zmin_se,
				reflector_depth, xmax_local, xmin_local, zmin_local, zmax_local, xmin_inj_local, xmax_inj_local, zmin_inj_local, zmax_inj_local, 
				xmax_se_local,  xmin_se_local, zmin_se_local):
	vp_max =np.max(SM_vp)
	vs_max =np.max(SM_vs)
	rho_max =np.max(SM_rho)
	vp_p = SM_vp
	vp_p[xmin_inj:xmax_inj, zmax-30:zmax-20]= vp_max
	vs_p = SM_vs
	vs_p[xmin_inj:xmax_inj, zmax-30:zmax-20]= vs_max
	rho_p = SM_rho
	rho_p[xmin_inj:xmax_inj,zmax-30:zmax-20]= rho_max
	
	shape   = (vp_p.shape[0],vp_p.shape[1])
		
	#Define index of the local domain, injection and extrapolation surface referencing the index in the full domain
	################
	#Local domain
	xmin=0
	xmax=500
	zmin=70
	zmax=150
	################
	#Injection
	xmin_inj= xmin
	xmax_inj= xmax
	zmin_inj= zmin+20
	zmax_inj= zmax
	################
	#Extrapolation surface (Se)
	xmin_se = 0.0
	xmax_se = 500
	zmin_se = zmin
	################
	#Reflector position in depth in the full domain 
	reflector_depth = zmin+30
	#################################################################################
	
	#################################################################################
	#Define index of the local domain, injection and extrapolation surface referencing the index in local domain
	#Local domain
	xmin_local=xmin-xmin
	xmax_local=xmax-xmin
	zmin_local=zmin-zmin
	zmax_local=zmax-zmin
	################
	#Injection
	xmin_inj_local = xmin_local
	xmax_inj_local = xmax_local
	zmin_inj_local = zmin_local+20
	zmax_inj_local = zmax_local
	#############
	#Extrapolation surface (Se)
	xmin_se_local= xmin_local
	xmax_se_local= xmax_local
	zmin_se_local= zmin_local
	#################
	#Reflector position in depth in local domain
	reflector_depth_local = zmin_local+30
	
	##################################################################################
	#Save full model 
	###################################################################################
	 
	outdict = dict()
	outdict['vp']     = vp_p
	outdict['vs']     = vs_p
	outdict['rho']    = rho_p
	outdict['spacing']= space
	outdict['shape']  = vp_p.shape
	outdict['nbl']    = nbl
	outdict['t0']     = t0
	outdict['tn']     = tn
	##Local domain index in the full domain
	outdict['xmin']=xmin
	outdict['xmax']=xmax
	outdict['zmin']=zmin
	outdict['zmax']=zmax
	##Injection full index
	outdict['xmin_inj']=xmin_inj
	outdict['xmax_inj']=xmax_inj
	outdict['zmin_inj']=zmin_inj
	outdict['zmax_inj']=zmax_inj
	##Extrapolation surface index
	outdict['xmin_se']   =xmin_se
	outdict['xmax_se']   =xmax_se
	outdict['zmin_se']   =zmin_se
	## Reflector position
	outdict['reflector_depth']   = reflector_depth
	
	## Local domain index
	outdict['xmax_local']=xmax_local
	outdict['xmin_local']=xmin_local
	outdict['zmin_local']=zmin_local
	outdict['zmax_local']=zmax_local
	##injection local index
	outdict['xmin_inj_local']=xmin_inj_local
	outdict['xmax_inj_local']=xmax_inj_local
	outdict['zmin_inj_local']=zmin_inj_local
	outdict['zmax_inj_local']=zmax_inj_local
	##extrapolation surface local index
	outdict['xmax_se_local']=xmax_se_local
	outdict['xmin_se_local']=xmin_se_local
	outdict['zmin_se_local']=zmin_se_local
	## Reflector position
	outdict['reflector_depth_local']   = reflector_depth_local
	
	spio.savemat('output/full_pert.mat', outdict)
	return vp_p


##########################################################################################

##load SEAM I 2D model
seam_vp =spio.loadmat('SEAM_I_2D_Model/SEAM_Vp.mat')
seam_den=spio.loadmat('SEAM_I_2D_Model/SEAM_Den.mat')
seam_vs =spio.loadmat('SEAM_I_2D_Model/SEAM_Vs.mat')

SEAM_vp_Elastic_N23900=seam_vp['SEAM_Vp']
SEAM_vs_Elastic_N23900=seam_vs['SEAM_Vs']
SEAM_den_Elastic_N23900=seam_den['SEAM_Den']

##Define grid spacing (x,z)
spacing = (20, 10) 
dx=20
dz=10

## cut the new model 
SM_vp  = (SEAM_vp_Elastic_N23900 [600:801, 1180:1681].T)/1000  ##[z,x]
SM_vs  = (SEAM_vs_Elastic_N23900 [600:801, 1180:1681].T)/1000
SM_rho = (SEAM_den_Elastic_N23900[600:801, 1180:1681].T)/1000

##QC plot local domain and extenden in the real coord 
# x_min=0.0
# z_min=0.0
# x_max=SM_vp.shape[0]*dx
# z_max=SM_vp.shape[1]*dz
# print(x_max,z_max)
# x_min_km = x_min/1000.0; x_max_km = x_max/1000.0; z_min_km = z_min/1000.0; z_max_km = z_max/1000.0; 
# plt.imshow(SM_vp[:].T, extent=[x_min_km, x_max_km,z_max_km, z_min_km])
# plt.title('subdomain SEAM Vp ')
# plt.xlabel('Distance x (km)')
# plt.ylabel('Depth (km)')
# plt.savefig('fig/Vp_subdomain.png')

space   = (20, 10)  # Grid spacing in m. The domain size is now 10km by 2km
origin  = (0., 0.)  # What is the location of the top left corner. This is necessary to define
nbl     = 200       # Number of points - absorption boundary condition
t0      = 0.0		# Initial time
tn      = 1300      # Final time
#################################################################################
#################################################################################

#Define index of the local domain, injection and extrapolation surface referencing the index in the full domain
################
#Local domain
xmin=0
xmax=500
zmin=70
zmax=150
################
#Injection
xmin_inj= xmin
xmax_inj= xmax
zmin_inj= zmin+20
zmax_inj= zmax
################
#Extrapolation surface (Se)
xmin_se = 0.0
xmax_se = 500
zmin_se = zmin
################
#Reflector position in depth in the full domain 
reflector_depth = zmin+30
#################################################################################

#################################################################################
#Define index of the local domain, injection and extrapolation surface referencing the index in local domain
#Local domain
xmin_local=xmin-xmin
xmax_local=xmax-xmin
zmin_local=zmin-zmin
zmax_local=zmax-zmin
################
#Injection
xmin_inj_local = xmin_local
xmax_inj_local = xmax_local
zmin_inj_local = zmin_local+20
zmax_inj_local = zmax_local
#############
#Extrapolation surface (Se)
xmin_se_local= xmin_local
xmax_se_local= xmax_local
zmin_se_local= zmin_local
#################
#Reflector position in depth in local domain
reflector_depth_local = zmin_local+30


vp, vs, rho, shape, spacing = full_model(nbl, t0, tn, space, SM_vp, SM_vs, SM_rho, xmin, xmax, zmin, zmax, xmin_inj, xmax_inj, zmin_inj, zmax_inj, xmin_se, xmax_se, zmin_se,
				reflector_depth, xmax_local, xmin_local, zmin_local, zmax_local, xmin_inj_local, xmax_inj_local, zmin_inj_local, zmax_inj_local, 
				xmax_se_local,  xmin_se_local, zmin_se_local)

# plt.figure()
# im=plt.imshow(vp.T)
# plt.savefig('fig/full.png')	


vp_p= full_model_pert(nbl, t0, tn, space, SM_vp, SM_vs, SM_rho, xmin, xmax, zmin, zmax, xmin_inj, xmax_inj, zmin_inj, zmax_inj, xmin_se, xmax_se, zmin_se,
				reflector_depth, xmax_local, xmin_local, zmin_local, zmax_local, xmin_inj_local, xmax_inj_local, zmin_inj_local, zmax_inj_local, 
				xmax_se_local,  xmin_se_local, zmin_se_local)

# plt.figure()
# im=plt.imshow(vp_p.T)
# plt.savefig('fig/full_p.png')

 
# so=4
# model = ModelElastic(vp=vp, vs=vs, b=1./rho, origin=origin, shape=shape, spacing=spacing, space_order=so, nbl=nbl)







	
	

