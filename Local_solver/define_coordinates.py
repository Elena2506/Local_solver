"""
    Local Solver

    This Python program define receivers, injection and extrapolation surface coordinates in full and local domain.

    The program is under the conditions terms in the file LICENSE.txt

    Authors: Ligia Elena Jaimes-Osorio, Alison Malcolm and Polina Zheglova, 2020.

    SEG version of this code available from:
    http://software.seg.org/2021/
"""
 
import numpy as np
import scipy.io as spio
from devito import *
from examples.seismic import *

###FULL DOMAIN 
################################################################################
##Load full model index 
#################################################################################
dict_fullmodel =spio.loadmat('output/full.mat')

vp      = dict_fullmodel['vp']
vs      = dict_fullmodel['vs']
rho     = dict_fullmodel['rho']
space   = dict_fullmodel['spacing']
nbl     = dict_fullmodel['nbl'][0][0]
t0      = dict_fullmodel['t0'][0][0]
tn      = dict_fullmodel['tn'][0][0]
origin  = (0., 0.) 
nbpml   = nbl
shape   = (vp.shape[0],vp.shape[1])
spacing = (space[0][0], space[0][1])
dx = space[0][0]
dz = space[0][1]


### load index of the local domain, injection and extrapolation surface referencing the full domain  ####
xmin=dict_fullmodel['xmin'][0][0]
xmax=dict_fullmodel['xmax'][0][0]
zmin=dict_fullmodel['zmin'][0][0]
zmax=dict_fullmodel['zmax'][0][0]

xmin_inj= dict_fullmodel['xmin_inj'][0][0]
xmax_inj= dict_fullmodel['xmax_inj'][0][0]
zmin_inj= dict_fullmodel['zmin_inj'][0][0]
zmax_inj= dict_fullmodel['zmax_inj'][0][0]

xmin_se = dict_fullmodel['xmin_se'][0][0]
xmax_se = dict_fullmodel['xmax_se'][0][0]
zmin_se = dict_fullmodel['zmin_se'][0][0]

## load index of the local domain, injection and extrapolation surface referencing the local domain  ####
xmax_local=dict_fullmodel['xmax_local'][0][0]
xmin_local=dict_fullmodel['xmin_local'][0][0]
zmin_local=dict_fullmodel['zmin_local'][0][0]
zmax_local=dict_fullmodel['zmax_local'][0][0]

xmin_inj_local=dict_fullmodel['xmin_inj_local'][0][0]
xmax_inj_local=dict_fullmodel['xmax_inj_local'][0][0]
zmin_inj_local=dict_fullmodel['zmin_inj_local'][0][0]
zmax_inj_local=dict_fullmodel['zmax_inj_local'][0][0]

xmax_se_local=dict_fullmodel['xmax_se_local'][0][0]
xmin_se_local=dict_fullmodel['xmin_se_local'][0][0]
zmin_se_local=dict_fullmodel['zmin_se_local'][0][0]

reflector_depth = dict_fullmodel['reflector_depth'][0][0]
######################################################################################

######################################################################################
###Define receiver coordinates 
######################################################################################

xmax_rec = (xmax*dx)                            #maximum receiver location in x
xmin_rec = (xmin*dx)                            #minimun receiver location in x
d_rec    = 10*dx                                #spacing between receivers
nrec     = int((xmax_rec-xmin_rec)/d_rec)       #number of receivers
z_rec    = dz                                   #receiver location in z

#####Save receiver coordinates
outdict              = dict()
outdict['xmax_rec']  = xmax_rec
outdict['xmin_rec']  = xmin_rec
outdict['nrec']      = nrec 
outdict['drec']      = d_rec
outdict['z_rec']     = z_rec
spio.savemat('output/coord_rec.mat', outdict)


######################################################################################
###Define surface extrapolation (Se) coordinates 
######################################################################################

x_max_se = (xmax_se*dx)                                     #maximum Se location in x
x_min_se = (xmin_se*dx)										#minimun Se location in x
d_se     = 2*dx												#spacing between Se points
n_se     = int((x_max_se-x_min_se)/d_se)					#number of Se points
z_se     = zmin_se*dz                                       #receiver location in z

### Define vector for coordinates
se_coordinates        = np.empty((n_se, 2)) 
se_coordinates[:,  0] = np.arange(x_min_se, x_max_se, d_se)
se_coordinates[:, -1] = z_se

#####Save Se points coordinates
outdict             = dict()
outdict['se_coord'] = se_coordinates
outdict['n_se']     = n_se
spio.savemat('output/coord_se.mat', outdict)


######################################################################################
###Define injection (S) coordinates 
######################################################################################

x_max_inj = xmax_inj*dx                                      #maximum S location in x
x_min_inj = xmin_inj*dx										 #minimun S location in x
dx_s=dx														 #spacing between S points in x
nx_inj    = int((x_max_inj-x_min_inj)/dx_s)                  #number of S points  in x

z_max_inj = zmax_inj*dz                                      #maximum S location in z
z_min_inj = zmin_inj*dz										 #minimun S location in z
dz_s=dz														 #spacing between S points in z
nz_inj=int((z_max_inj-z_min_inj)/dz_s)                       #number of S points  in z

### Define vectors for coordinates
###top
inj_coordinates_t        = np.empty((nx_inj, 2)) 
inj_coordinates_t[:, 0]  = np.arange(x_min_inj, x_max_inj, dx_s) 
inj_coordinates_t[:, -1] = z_min_inj

###bot
inj_coordinates_b        = np.empty((nx_inj, 2))  
inj_coordinates_b[:, 0]  = np.arange(x_min_inj, x_max_inj, dx_s)  
inj_coordinates_b[:, -1] = z_max_inj

##left
inj_coordinates_l = np.empty((nz_inj, 2)) 
inj_coordinates_l[:, 0] = x_min_inj
inj_coordinates_l[:, -1]= np.arange(z_min_inj, z_max_inj, dz_s)
 
##right
inj_coordinates_r = np.empty((nz_inj, 2))  
inj_coordinates_r[:, 0] = x_max_inj
inj_coordinates_r[:, -1]= np.arange(z_min_inj, z_max_inj, dz_s)

### concatene 
inj_t = np.concatenate([inj_coordinates_t, inj_coordinates_l, inj_coordinates_b, inj_coordinates_r])

#####Save S points coordinates
outdict                     = dict()
outdict['inj_coordinates_t']= inj_coordinates_t
outdict['inj_coordinates_l']= inj_coordinates_l
outdict['inj_coordinates_b']= inj_coordinates_b
outdict['inj_coordinates_r']= inj_coordinates_r
outdict['inj_t']            = inj_t
outdict['nx_inj']  			=nx_inj
outdict['nz_inj']			=nz_inj
spio.savemat('output/coord_inj.mat', outdict)


######################################################################################
###Define source coordinates 
######################################################################################

xsrc    = [xmax*dx/2]          # x source position
nsrc    = len(xsrc)            # number of sources 
zsrc    = dz                   # z source position 
dsrc   = xmax*dx/2			   # spacing between src points

src_coordinates = np.empty((nsrc, 2))
src_coordinates[:,  0] = xsrc
src_coordinates[:, -1] = zsrc

######Save source coordinates
outdict          = dict()
outdict['xsrc']  = xsrc
outdict['zsrc']  = zsrc
outdict['nsrc']  = nsrc 
outdict['dsrc']  = dsrc
outdict['src_coord']=src_coordinates
spio.savemat('output/coord_src.mat', outdict)
################################################################################


###LOCAL DOMAIN 
################################################################################
##Load local model index 
#################################################################################
dict_model =spio.loadmat('output/local.mat')

vp      = dict_model['vp']
vs      = dict_model['vs']
rho     = dict_model['rho']
space   = dict_model['spacing']
nbl     = dict_model['nbl'][0][0]
t0      = dict_model['t0'][0][0]
tn      = dict_model['tn'][0][0]
origin  = (0., 0.) 
nbpml   = nbl
shape   = (vp.shape[0],vp.shape[1])
spacing = (space[0][0], space[0][1])
dx = space[0][0]
dz = space[0][1]


### load index of the local domain, injection and extrapolation surface referencing the full domain  ####
xmin=dict_model['xmin'][0][0]
xmax=dict_model['xmax'][0][0]
zmin=dict_model['zmin'][0][0]
zmax=dict_model['zmax'][0][0]

xmin_inj= dict_model['xmin_inj'][0][0]
xmax_inj= dict_model['xmax_inj'][0][0]
zmin_inj= dict_model['zmin_inj'][0][0]

xmin_se = dict_model['xmin_se'][0][0]
xmax_se = dict_model['xmax_se'][0][0]
zmin_se = dict_model['zmin_se'][0][0]


## load index of the local domain, injection and extrapolation surface referencing the local domain  ####
xmax_local_l=dict_model['xmax_local'][0][0]
xmin_local_l=dict_model['xmin_local'][0][0]
zmin_local_l=dict_model['zmin_local'][0][0]
zmax_local_l=dict_model['zmax_local'][0][0]

xmin_inj_l=dict_model['xmin_inj_local'][0][0]
xmax_inj_l=dict_model['xmax_inj_local'][0][0]
zmin_inj_l=dict_model['zmin_inj_local'][0][0]
zmax_inj_l=dict_model['zmax_inj_local'][0][0]

xmax_se_l=dict_model['xmax_se_local'][0][0]
xmin_se_l=dict_model['xmin_se_local'][0][0]
zmin_se_l=dict_model['zmin_se_local'][0][0]
######################################################################################

######################################################################################
###Define injection (S) coordinates local
######################################################################################


x_max_inj = xmax_inj_l*dx                                   #maximum S location in x               
x_min_inj = xmin_inj_l*dx									#minimun S location in x
dx_s      = dx												#spacing between S points in x
nx_inj    = int((x_max_inj-x_min_inj)/dx_s)                 #number of S points in x

z_max_inj = zmax_inj_l*dz                                   #maximum S location in z
z_min_inj = zmin_inj_l*dz									#minimun S location in z
dz_s      = dz
nz_inj    = int((z_max_inj-z_min_inj)/dz_s)                 #number of S points  in z

### Define vectors for coordinates
###top
inj_coordinates_t        = np.empty((nx_inj, 2)) 
inj_coordinates_t[:, 0]  = np.arange(x_min_inj, x_max_inj, dx_s) 
inj_coordinates_t[:, -1] = z_min_inj

###bot
inj_coordinates_b        = np.empty((nx_inj, 2))  
inj_coordinates_b[:, 0]  = np.arange(x_min_inj, x_max_inj, dx_s)  
inj_coordinates_b[:, -1] = z_max_inj

##left
inj_coordinates_l = np.empty((nz_inj, 2)) 
inj_coordinates_l[:, 0] = x_min_inj
inj_coordinates_l[:, -1]= np.arange(z_min_inj, z_max_inj, dz_s)
 
##right
inj_coordinates_r = np.empty((nz_inj, 2))  
inj_coordinates_r[:, 0] = x_max_inj
inj_coordinates_r[:, -1]= np.arange(z_min_inj, z_max_inj, dz_s)

### concatene 
inj_t = np.concatenate([inj_coordinates_t, inj_coordinates_l, inj_coordinates_b, inj_coordinates_r])

#####Save S points coordinates local 
outdict       				= dict()
outdict['inj_coordinates_t']= inj_coordinates_t
outdict['inj_coordinates_l']= inj_coordinates_l
outdict['inj_coordinates_b']= inj_coordinates_b
outdict['inj_coordinates_r']= inj_coordinates_r
outdict['inj_t']			= inj_t
outdict['nx_inj']			= nx_inj
outdict['nz_inj']			= nz_inj
spio.savemat('output/coord_inj_local.mat', outdict)


######################################################################################
###Define surface extrapolation (Se) coordinates 
######################################################################################

x_max_se = (xmax_se_l*dx)                                   #maximum Se location in x
x_min_se = (xmin_se_l*dx)									#minimun Se location in x
z_min_se = (zmin_se_l*dz)                                   #receiver location in z
d_se     = 2*dx												#spacing between Se points
n_se=int((x_max_se-x_min_se)/d_se)							#number of Se points

### Define vector for coordinates
se_coordinates        = np.empty((n_se, 2)) 
se_coordinates[:,  0] = np.arange(x_min_se, x_max_se, d_se)
se_coordinates[:, -1] = z_min_se

#####Save Se points coordinates
outdict            = dict()
outdict['se_coord']= se_coordinates
outdict['n_se']    = n_se
spio.savemat('output/coord_se_local.mat', outdict)
