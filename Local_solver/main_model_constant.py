"""
    Local Solver

    A Python main code to run Example 1 constant elastic velocity model. 
    
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
###EXAMPLE 1. Constant velocity model 
################################################################################

## Ram memory require for all the sequence ~ 30Gbytes. Time ~ 30 minutes 
## To run this main code you need:

## 1. Create a model
## Run file create_model_constant.py


## 2. Define coordinates
## Run file define_coordinates.py

## After steps 1 and 2 you are able to run this main code. 
## 3. Define input parameter
## Define space-order=so

so=4
## 4. Load source positions
dict_src =spio.loadmat('output/coord_src.mat') 
src_coordinate = dict_src['src_coord']

## 5. Forward modelling and local injection simulation on initial model
full_simulation (src_coordinate, so, model_dirt='output/full.mat', mi='nonp')
local_simulation(src_coordinate, so, model_dirt='output/local.mat',mi='nonp')

## 6. Plot model and full and injection comparison 
geom='output/geometry_nonp_%s.00.mat'%int(src_coordinate[0][0])
plot_model(model='output/full.mat', geom=geom)
plot_injection(geom, v_local='output/v1_local_nonp.mat', model = 'output/full.mat', mo='nonp')
plot_simple_traces (arr_full='output/waveform_se_full_src_nonp_2500.00.mat', arr_local='output/waveform_se_local_src_nonp_2500.00.mat', xp=20, time_i=380, time_f=440 )

################################################################################
