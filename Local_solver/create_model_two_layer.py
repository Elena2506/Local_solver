"""
    Local Solver

    A Python program to:
    - Compute two layer velocity model.
    - Define Full and Local domain.  

    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""


import numpy as np
import scipy.io as spio
from devito import *
from examples.seismic import *

#################################################################################
##DEFINE full MODEL #########
#################################################################################
# Define a physical size 

shape   = (501, 201)  # Number of grid point (nx, nz)
spacing = (10., 10.)  # Grid spacing in m. The domain size is now 1km by 1km
origin  = (0., 0.)    # What is the location of the top left corner. This is necessary to define

nbl     = 200         # Number of points - absorption boundary condition
t0      = 0.0		   # Initial time
tn      = 1700         # Final time

# Define a velocity profile. The velocity is in km/s
vp = np.empty(shape, dtype=np.float32)
vp [:]               = 2.0
vp [:, 100:]         = 4.0
vp[500:501,200:201] = 5.3
vs = np.empty(shape, dtype=np.float32)
vs [:]               = 0.88
vs [:, 100:]         = 1.54
vs[500:501,200:201]  = 2.35
rho = np.empty(shape, dtype=np.float32)
rho[:]               = 2.0
rho[:, 100:]         = 2.3
rho[500:501,200:201] = 2.86
##################################################################################

#Define index of the local domain, injection and extrapolation surface referencing the index in the full domain
################
#Local domain
xmin=0
xmax=500
zmin=140
zmax=190
################
#Injection
xmin_inj= xmin
xmax_inj= xmax
zmin_inj= zmin+10
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
zmin_inj_local = zmin_local+10
zmax_inj_local = zmax_local
#############
#Extrapolation surface (se)
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
outdict['vp']     = vp
outdict['vs']     = vs
outdict['rho']    = rho
outdict['spacing']= spacing
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
##Extrapolation surface  index
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

spio.savemat('output/full.mat', outdict)

#################################################################################
##Define full model perturbed #########
#################################################################################

# Define a velocity profile. The velocity is in km/s
vp_p = np.empty(shape, dtype=np.float32)
vp_p [:]                                 = 2.0
vp_p [:,100:]                            = 4.0
vp_p[xmin_inj:xmax_inj,zmin+30:zmax-10]  = 5.3
vs_p = np.empty(shape, dtype=np.float32)
vs_p [:]                                 = 0.88
vs_p [:, 100:]                           = 1.54
vs_p[xmin_inj:xmax_inj,zmin+30:zmax-10]  = 2.35
rho_p = np.empty(shape, dtype=np.float32)
rho_p[:]                                 = 2.0
rho_p[:, 100:]                           = 2.3
rho_p[xmin_inj:xmax_inj,zmin+30:zmax-10] = 2.86

##################################################################################
#Save full model perturbed
###################################################################################
 
outdict = dict()
outdict['vp']      = vp_p
outdict['vs']      = vs_p
outdict['rho']     = rho_p
outdict['spacing'] = spacing
outdict['shape']   = vp.shape
outdict['nbl']     = nbl
outdict['t0']      = t0
outdict['tn']      = tn

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
##Extrapolation surface  index
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


#################################################################################
##Define local domain #########
#################################################################################

# Define a physical size
shape_lnonp =  (xmax_local+1, zmax_local+1)   # Number of grid point (nx, nz)
spacing_lnonp = (10., 10.)                    # Grid spacing in m. The domain size is now 1km by 1km
origin = (0., 0.)  							  # What is the location of the top left corner. This is necessary to define

# Define a velocity profile. The velocity is in km/s
vp_lnonp = np.empty(shape_lnonp, dtype=np.float32)
vp_lnonp[:]                                               = 4.0
vp_lnonp[xmax_local:xmax_local+1,zmax_local:zmax_local+1] = 5.3
vs_lnonp = np.empty(shape_lnonp, dtype=np.float32)
vs_lnonp[:]                                               = 1.54
vs_lnonp[xmax_local:xmax_local+1,zmax_local:zmax_local+1] = 2.35
rho_lnonp = np.empty(shape_lnonp, dtype=np.float32)
rho_lnonp[:]                                               = 2.3
rho_lnonp[xmax_local:xmax_local+1,zmax_local:zmax_local+1] = 2.86

##################################################################################
#Save local model
###################################################################################
 
outdict = dict()
outdict['vp']      = vp_lnonp
outdict['vs']      = vs_lnonp
outdict['rho']     = rho_lnonp
outdict['spacing'] = spacing_lnonp
outdict['shape']   = vp_lnonp.shape
outdict['nbl']     = nbl
outdict['t0']      = t0
outdict['tn']      = tn

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
##Extrapolation surface  index
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

spio.savemat('output/local.mat', outdict)

#################################################################################
##Define local domain perturbed #########
#################################################################################

# Define a physical size
shape_lp = (xmax_local+1, zmax_local+1)  # Number of grid point (nx, nz)
spacing_lp = (10., 10.)                  # Grid spacing in m. The domain size is now 1km by 1km
origin = (0., 0.)                        # What is the location of the top left corner. This is necessary to define

# Define a velocity profile. The velocity is in km/s
vp_lp    = np.empty(shape_lp, dtype=np.float32)
vp_lp[:]                                                     = 4.0
vp_lp[xmin_local:xmax_local+1,zmin_local+30:zmax_local-10]   = 5.3
vs_lp = np.empty(shape_lp, dtype=np.float32)
vs_lp[:]                                                     = 1.54
vs_lp[xmin_local:xmax_local+1,zmin_local+30:zmax_local-10]   = 2.35
rho_lp = np.empty(shape_lp, dtype=np.float32)
rho_lp[:]                                                    = 2.3
rho_lp[xmin_local:xmax_local+1,zmin_local+30:zmax_local-10]  = 2.86

#################################################################################
# SAVE DICT MODEL LOCAL PERTU
#################################################################################


outdict = dict()
outdict['vp']        = vp_lp
outdict['vs']        = vs_lp
outdict['rho']       = rho_lp
outdict['spacing']   = spacing_lp
outdict['shape']     = vp_lp.shape
outdict['nbl']       = nbl
outdict['t0']        = t0
outdict['tn']        = tn

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
##Extrapolation surface  index
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

spio.savemat('output/local_pert.mat', outdict)



