"""
    Local Solver

    A Python program to plot velocity model seam
    
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


###load the model 
dict_submodel =spio.loadmat('output/full_pert.mat')

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
dict_src =spio.loadmat('output/coord_src.mat') 
src_coordinate = dict_src['src_coord']
geom='output/geometry_nonp_%s.00.mat'%int(src_coordinate[0][0])
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
plt.subplot(311)
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

plt.colorbar(label='V(km/s)')
plt.title('P-velocity full' )
# plt.xlabel('Horizontal coordinate (km)')
plt.ylabel('Depth (km)')  

plt.subplot(312)
arr=rho.T
im = plt.imshow(arr, extent=[x_min_km, x_max_km,z_max_km, z_min_km])
plt.colorbar(label='Density(kg/m3)')
plt.title('Density full' )
# plt.xlabel('Horizontal coordinate (km)')
plt.ylabel('Depth (km)')


plt.subplot(313)
arr=vs.T
im = plt.imshow(arr, extent=[x_min_km, x_max_km,z_max_km, z_min_km])
plt.colorbar(label='V(km/s)')
plt.title('S-velocity full' )
plt.xlabel('Horizontal coordinate (km)')
plt.ylabel('Depth (km)')
plt.tight_layout()
plt.savefig('fig/model_full.png')


















