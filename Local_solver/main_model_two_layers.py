"""
    Local Solver

    A Python main code to run Example 1 Constant elastic velocity model. 
    
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
from full_simulation import *
from local_simulation import *
from plot import *
from compute_green_function import *
from convolution import convolution 



################################################################################
## main code		
################################################################################
###EXAMPLE 2. Two layer velocity model 
################################################################################
 
## Ram memory require for all the sequence ~ 30Gbytes. Time ~ 30 minutes 
## To run this main code you need:

## 1. Create a model
## Run file create_model_two_layer.py


## 2. Define coordinates
## Run file define_coordinates.py

## After steps 1 and 2 you are able to run this main code. 
## 3. Define input parameter
## Define space-order=so
so=4
## 4. Load source positions. For now we use only one shot
dict_src =spio.loadmat('output/coord_src.mat') 
src_coordinate = dict_src['src_coord']

## 5. Forward modelling and local injection simulation on initial model
full_simulation (src_coordinate, so, model_dirt='output/full.mat', mi='nonp')
local_simulation(src_coordinate, so, model_dirt='output/local.mat',mi='nonp')
 
## 6. Plot model and full and injection comparison 
geom='output/geometry_nonp_%s.00.mat'%int(src_coordinate[0][0])
plot_injection(geom, v_local='output/v1_local_nonp.mat', model = 'output/full.mat', mo='nonp')

## 7. Forward modelling and local injection simulation on perturbed model 
full_simulation (src_coordinate, so, model_dirt='output/full_pert.mat',  mi='pert')
local_simulation(src_coordinate, so, model_dirt='output/local_pert.mat', mi='pert')

## 8. Plot model and full and injection comparison 
plot_model(model='output/full_pert.mat', geom=geom)
plot_injection( geom='output/geometry_pert_%s.00.mat'%int(src_coordinate[0][0]), v_local='output/v1_local_pert.mat', model = 'output/full_pert.mat', mo='pert')

## 9. Compute Green functions in each direction 
Green_function(model='output/full_pert.mat', geom=geom, so=so, direction=0, save_dir='x')
Green_function(model='output/full_pert.mat', geom=geom, so=so, direction=1, save_dir='z')

## 10. Convolution 
convolution(geom, wave='output/waveform_se_local_src_pert_%s.00.mat'%int(src_coordinate[0][0]))

## 11. Plot shotgather 
plot_shotgather(shotgather_extra='output/shotgather_extrapolated.mat', shotgather_pert ='output/shotgather_src_pert_%s.00.mat'%int(src_coordinate[0][0]), geom=geom)

## 12. Plot traces
plot_traces(shotgather_extra='output/shotgather_extrapolated.mat', geom=geom, shotgather_wf= 'output/shotgather_wf.mat', xp=20, time_i=1240, time_f=1340) 

# 13. Amplitude and Phase extraction and comparison
#Run file amplitude_phase.py  
################################################################################


