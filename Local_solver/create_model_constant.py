"""
    Local Solver

    A Python program to:
    Create constant elastic velocity model.
    Define Full and Local domain.  

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
shape   = (501, 201)   # Number of grid point (nx, nz)
spacing = (10., 10.)   # Grid spacing in m. The domain size is now 1km by 1km
origin  = (0., 0.)     # What is the location of the top left corner. This is necessary to define

nbl     = 100          # Number of points - absorption boundary condition
t0      = 0.0		   # Initial time
tn      = 2000         # Final time

# Define velocity profile. The velocity is in km/s
vp = np.empty(shape, dtype=np.float32)
vp [:]= 2.0
vs = np.empty(shape, dtype=np.float32)
vs [:]= 0.88
rho = np.empty(shape, dtype=np.float32)
rho[:] = 2.0

##################################################################################
#Define index of the local domain, injection and extrapolation surface referencing the index in the full domain
#Local domain
xmin=0
xmax=500
zmin=100
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
zmin_se = zmin+15
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
zmin_se_local= zmin_local+15
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

##################################################################################
##################################################################################


##################################################################################
#################################################################################
##Define local domain #########
#################################################################################

# Define a physical size
shape_lnonp =  (xmax_local+1, zmax_local+1)   # Number of grid point (nx, nz)
spacing_lnonp = spacing                       # Grid spacing in m. 
origin = (0., 0.)                             # Top left corner in the model.


# Define a velocity profile. The velocity is in km/s
vp_lnonp     = np.empty(shape_lnonp, dtype=np.float32)
vp_lnonp[:]  = 2.0
vs_lnonp     = np.empty(shape_lnonp, dtype=np.float32)
vs_lnonp[:]  = 0.88
rho_lnonp    = np.empty(shape_lnonp, dtype=np.float32)
rho_lnonp[:] = 2.0


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


