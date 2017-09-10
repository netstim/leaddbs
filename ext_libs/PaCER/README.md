## Welcome to PaCER - Precise and Convenient Electrode Reconstruction for DBS

Please note that PaCER is a research tool **NOT** intended for clinical use.   
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**GNU Affero General Public License** for more details.

Copyright (C) 2016-2017  Andreas Husch, 
					Centre Hospitalier de Luxembourg, National Department of Neurosurgery and
				    University of Luxembourg, Luxembourg Centre for Systems Biomedicine
![Image of a PaCER electrode reconstruction at two different time points of resolving brain shift.](docs/PaCER.png)
### Background 
The PaCER Toolbox is a MATLAB implementation of a robust method to fully automatically reconstruct deep brain stimulation trajectories from post operative CT imaging data. PaCER is able to fully preserving electrode bending (e.g. caused by brainshift). Further is able to detect individual contacts on high-resolution data. 
The PaCER toolbox is provided with means to easily visualize electrodes as well as imaging data within the MATLAB environment
. 
### Requirements
The requirements to use PaCER are 
 *  **MATLAB**  
 *  a **post-operative CT image** in **nifti** file format. 
 
A **CT slice-thickness <= 1 mm** is recommend, however PaCER will
work on lower resolution data too by falling back to a less sophisticated contact detection 
method (yielding lower accuracy). Nifti input files are supported in compressed form (.nii.gz) as
well as non-compressed (.nii).

### Getting Started
The easiest way to learn about PaCER is to run the example files. We recommend to add the
PaCER directory and all its subdirectories to your MATLAB path first. This can be 
archived by running the file SETUP_PACER.m in MATLAB (once). The examples include a call
to SETUP_PACER.

#### Example Dataset
The following examples require only a post op CT dataset - they should work out of the box for most CT scan protocols as long as the slice thickness is not toooo bad :-) Easy conversion from DICOM to NIFTI is possible dcm2nii which is included in [MRIcron](https://www.nitrc.org/projects/mricron/).
Advanced example demonstrating further use-cases (e.g. visualisation of segmentations and atlas data, simple volume of tissue activated model etc.) can be found in the examples/advanced/ directory. However these examples require appropriate co-registered image modalities (e.g. atlases, segmentation). We are in the process to provide a full example dataset in the future.

#### The Examples
 * **EXAMPLE_1.m** - Basic PaCER call and electrode plot. Start here!
     * **EXAMPLE_1_1.m** - Continues EXAMPLE_1 by adding an **MPR view** of the CT image and demonstrating some **plot customisations**

### Questions
Feel free to open an issue at [https://github.com/adhusch/PaCER](https://github.com/adhusch/PaCER) or drop a note to mail (at) andreashusch.de

### Literature
A paper describing the algorithmic details of PaCER as well as reporting detailled accuracy analysis results on phantom and simulated data is currently in revision. It will be linked here after publication.

### Acknowledgement
This work was made possible by a Aide à la Formation Recherche grant (AFR) to Andreas Husch by the Luxembourg National Research (FNR).


PaCER is packaged with some free external software libraries for convenience. Please see the "toolboxes" folder and the respective LICENSE files for details.
We feel grateful to the authors of this toolboxes and scripts:
 * [Tools for NIfTI and ANALYZE image](https://de.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image) by Jimmy Shen
 * [RGB triple of color name, version 2](https://de.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name--version-2) by Kristjan Jonasson
 * [GUI Layout Toolbox](https://de.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox) by David Sampson and Ben Tordoff
 * [in_polyhedron](https://de.mathworks.com/matlabcentral/fileexchange/48041-in-polyhedron) by Jaroslaw Tuszynski
 * [Cylinder Between 2 Points](https://de.mathworks.com/matlabcentral/fileexchange/5468-cylinder-between-2-points) by Per Sundqvist
 * MPR View, by Florian Bernard

