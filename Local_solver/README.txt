README

SEG code 2020-0001, “Building a fast elastic solver”

By Ligia Elena Jaimes-Osorio  & Alison Malcolm and Polina Zheglova (Memorial University of Newfoundland)

This repository contains the manuscript and a collection of Python codes, synthetic data example and data file for a paper about the open-source software package Local-solver. The example data reproduces the results and figures shown in the SEG publication. The Local-solver software is written in Python 3.6.10 and programming languages. To run the program, the Devito, numpy, scipy and matplotlib library is required. (Devito: https://www.devitoproject.org/devito/). We recommend using Docker installation as is provided in (https://www.devitoproject.org/devito/download.html) After Docker and Devito are installed  this repository should be copied into the docker container, then a bash shell with Devito should be started using: "docker-compose run devito /bin/bash" before the code can be executed.



1-Abstract
2-Content
3-Prerequisites
4.Running the files
5-License
6-Disclamer

==========================================================

1- Abstract 
—————————————————

The recovery of elastic properties from seismic data often requires iterative the use of seismic modeling. The finite-difference method is an essential component in seismic modeling and usually the most computationally expensive step in methodologies such as inversion or reverse time migration. We introduce an implementation  of an elastic local solver that allows us to propagate the elastic wavefield within a subvolume after local alteration of the model. We implement the elastic local solver as an add-on to an existing solver, so that minimal changes are required to turn any finite-difference solver into a local solver. We demonstrate the capability of the elastic local solver to reconstruct the wavefield and then we extrapolate the wavefield to the receiver positions to compare with observed data. 


2- Content
—————————————————
Python codes

-General python modules containing the functions to calculate the shotgather, injection, wavefields, Green functions, and convolution:
	-full_simulation.py
	-local_simulation.py 
	-convolution.py
	-compute_green_function.py
	-amplitude_phase.py 
	-calc_phase.py

- Python scripts to run the different models:
	-main_model_constant.py
	-main_model_seam.py
	-main_model_two_layers.py

Python scripts to generate the synthetic local and full models from SEAM_I_2D_Model:
-create_model_seam_full.py
-create_model_seam_local.py
	    
Python script to generate the synthetic local and full constant models:
-create_model_constant.py
		
Python script to generate the synthetic local and full two layer models:
-create_model_two_layer.py
		
Python script to generate the synthetic local and full two layer coordinates:
-define_coordinates.py

Python general modules that contains functions to plot: model, injection, shotgather, and traces:			
-plot.py
-plot_model.py		 

Data files: 
Synthetic SEAM data. This data can be download: https://seg.org/News-Resources/Research-Data/Open-Data 
-SEAM_I_2D_Model

3- Prerequisites
—————————————————
The Python software Local solver - requires the library “Devito”, and some Python packages.

The easier way to get Devito is following the instruction 

https://www.devitoproject.org/devito/download.html

After installed Devito, the software “Local solver” can run in terminal. 

The “Local solver” programs were tested Mac platforms using Python 3.6.10.

Ram memory require for all the sequence ~ 30Gbytes. Time ~ 30 minutes 

4- Running the files
—————————————————
After the libraries are installed, download the folder “Local_solver”, on your machine. You can download a zip archive of this
repository at: software.seg.org(http://software.seg.org). Following, be sure that all files are in the same directory. 

Open the command prompt window (terminal) and run the command:

1. Constant model:
- python create_model_constant.py
- python define_coordinates.py
- python main_model_constan.py

2. Two layer model:  
- python create_model_two_layer.py
- python define_coordinates.py
- Two layer model:  python main_model_two_layers.py
- Two layer model:  python amplitude_phase.py

3. SEAM model:

- python create_model_seam_full.py
- python create_model_seam_local.py
- python define_coordinates.py
- python main_model_seam.py 
- python model_full.png

WARNING: The output files overwrite when you change the model. Make sure you move the outputs before start to run a new model. 

5- Output files
—————————————————
1. File create_model_constant.py: 
-full.mat
-local.mat

2. File create_model_two_layers.py: 
-full.mat
-local.mat
-full_pert.mat
-local_pert.mat	

3. File python create_model_seam_full.py: 
-full.mat
-full_pert.mat

4.File python create_model_seam_local.py:
-local.mat
-local_pert.mat	

5.File define_coordinates.py:
-coord_rec.mat
-coord_se.mat
-coord_inj.mat
-coord_inj_local.mat
-coord_se_local.mat

6.File full_simulation.py:
## mi= pert or nonp, isrc=source x coord. 
-shotgather_src_%s_%.2f.mat'%(mi, isrc)
-waveform_se_full_src_%s_%.2f.mat'%(mi, isrc)
-point_sources_src_%s_%.2f.mat'%(mi, isrc) #Point sources equation (2)
-geometry_%s_%.2f.mat'%(mi, isrc) 

7.File local_simulation.py
## mi= pert or nonp, isrc=source x coord. 
-waveform_se_local_src_%s_%.2f.mat'%(mi,i[0]) 
-v1_local_%s.mat'%mi

8.File convolution.py 
-shotgather_extrapolated.mat

9. compute_green_functions.py
## isorc=source x coord. 
-Component_source_greens_se%.2f'%isorc


10. main_model_constant.py:
-figure:model-vp.png 
-figure:injection_non_src_2500.00.png (figure 3 in the paper)
-figure:trace=%.2f.png'%(xp) (figure 4 in the paper) xp =receiver index position 

11. main_model_two_layers.py:

-figure: injection_nonp_src_2500.00.png
-figure: injection_pert_src_2500.00.png
-figure: model-vp.png (figure 6 in the paper)
-figure: shotgather_vx.png (figure 7 in paper)
-figure: shotgather_vz.png
-figure: vx_short_trace=20.00.png (figure 8 in paper)
-figure: vz_short_trace=20.00.png
-figure: vx_trace=20.00.png
-figure: vz_trace=20.00.png
-figure: vx_out_integral.png
-figure: vz_out_integral.png

12. amplitude_and_phase.py
-phases.mat
-figure: windowed_shotgather_true_vx.png
-figure: windowed_shotgather_extr_vx.png
-figure: amplitude.png (figure 9 in the paper)
-figure: Phase.png (figure 10 in the paper)
-figure: QC_shot_timeindex.png

13. main_model_seam.py
-figure: injection_nonp_src_2500.00.png
-figure: injection_pert_src_2500.00.png
-figure: model-vp.png 
-figure: shotgather_vx.png (figure 12 in the paper)
-figure: shotgather_vz.png
-figure: vx_short_trace=20.00.png (figure 13 in the paper)
-figure: vz_short_trace=20.00.png
-figure: vx_trace=20.00.png
-figure: vz_trace=20.00.png

14. plot_model.py
-figure: model_full.png (figure 11 in the paper)

6- License
—————————————————
The following legal note is restricted solely to the content of the named files. It cannot
overrule licenses from the Python standard distribution modules, which are imported and
used therein.

The “Local solver” folder are distributed under the following license
agreement:

Copyright (c) 2020, Ligia Elena Jaimes-Osorio & Alison Malcolm & Polina Zheglova (Memorial University of Newfoundland)


All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials provided 
with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
OF THE POSSIBILITY OF SUCH DAMAGE.


7- Disclaimer 
—————————————————
If you obtained this set of codes from the SEG (downloaded from software.seg.org or
otherwise), you must also agree to the following disclaimer:

http://software.seg.org/disclaimer2.txt