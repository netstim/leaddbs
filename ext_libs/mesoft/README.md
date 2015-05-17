# README #

### Idea

This projects implements the algorithm described in 

**MesoFT: Unifying Diffusion Modelling and Fiber Tracking**

M. Reisert, V.G. Kiselev, B. Dihtal, E. Kellner, and D.S. Novikov, 

MICCAI 2014

### What do you need? ###

* A running MATLAB version
* A c-compiler to compile c-files within MATLAB with mex
* If you want to use the parallel implementation, you need openMP
* The MATLAB Nifti-toolbox by Jimmy Shen to enable Nifti im/export, you can find it here [(here)](http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)

### How do I get set up? ###

* Get the sources and unpack it somewhere
* Within MATLAB change into the project folder
* Type *>> init_mesoft*. This sets all the paths and compiles the basic mex-files
* That's it, type *>> mesoGT_tool* to get started
* Use *>> help mesoGT_tool* for more information

### What kind of data can be fed in ###

Currently there are four ways.

1. You can pass simple matlab arrays to the tracker  via *mesoGT_tool('setData',data,tensor,mask,vox,edges,name)*
2. You can read .mat files in format provided by [DTI-Fiber Tools](http://www.uniklinik-freiburg.de/mr-en/research-groups/diffperf/fibertools.html) with *>> mesoGT_tool('loadData','mat',....)*
3. You can read Niftis in FSL-type format with *>> mesoGT_tool('loadData','nii',....)*
4. You can directly read HCP-data

### A simple phantom example ###

* There is a .mat file *example_2shell_scheme.mat* which contains a simple 2-shell q-space scheme, from this you can synthesize a numerical phantom to test the algorithm. 
* Load the scheme via *>> load example_2shell_scheme.mat* 
* Generate the phantom via *mesoGT_tool('genPhantom',ten,0.05)* at an SNR level of 20 = 1/0.05
* Press the 'Start Tracking' Button.
* Check for 'More Statistics' to get some nice visulization of tracking process

### Misc ###

* Note, that for every q-shell scheme the m-file mesoGT_tool creates some custom c-code which is compiled on-the-fly and stored in some temporary folder.
* To get an overview over the parameters check GTdefaults.m, there you can basically all parameters you can modify. You can also choose there which parameters should be shown in GUI. You can also change the parameters via mesoGT_tool('setparam',....)
* The tracts are currently saved in some custom .mat format which is readable by [DTI-Fiber Tools](http://www.uniklinik-freiburg.de/mr-en/research-groups/diffperf/fibertools.html) 
* You can save parametric maps of the diffusion parameters via the *Save PM*-button as Niftis
* and ....

### File descriptions ###