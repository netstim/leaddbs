## Welcome to PaCER - Precise and Convenient Electrode Reconstruction for DBS

Please note that PaCER is a research tool **NOT** intended for clinical use.   
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

Copyright (C) 2017  Andreas Husch, 
					Centre Hospitalier de Luxembourg, National Department of Neurosurgery and
				    University of Luxembourg, Luxembourg Centre for Systems Biomedicine

### Requirements
The requirements to use PaCER are 
 *  a working **MATLAB installation**  
 *  **post-operative CT image** in **nifti** file format. 
 
A CT slice-thickness <= 1 mm is recommend, however PaCER will
work on lower resolution data too by falling back to a less sophisticated contact detection 
method (yielding lower accuracy). Nifti input files are supported in compressed form (.nii.gz) as
well as non-compressed (.nii).

### Getting Started
The easiest way to learn about PaCER is to run the example files. We recommend to add the
PaCER directory and all its subdirectories to your MATLAB path first. This can be 
archived by running the file SETUP_PACER.m in MATLAB (once). The examples include a call
to SETUP_PACER.

#### The Examples
 * **EXAMPLE_1.m** - Basic PaCER call and electrode plot. Start here!
    * **EXAMPLE_1_1.m** - Continues EXAMPLE_1 by adding an **MPR view** of the CT image and demonstrating some **plot customisations**
    * **EXAMPLE_1_2.m** - Continues EXAMPLE_1 by adding **transformation of the electrode object to native T1 space** and loading MPR view of the native T1 image. This transformation is achieved by applying a FSL FLIRT transformation (from CT to T1 space) to the electrode model using the methods provided by the PaCER electrode objects.
 * **EXAMPLE_2.m** - Demonstrates the electrodes reconstruction from **longitudinal datasets** (immediate post-op, later post-op) as well as **visualisation of co-registered MRI images** (T1) as well as **visualisation of STN segmentations**. (cf. Use-Case A, Use-Case C).
 	* **EXAMPLE_2_1.m** - Continues EXAMPLE_2 by adding a **simplified VTA model** (Mädler/Coenen) (cf. Use-Case C)
 * **EXAMPLE_3.m** - Demonstrates plan/outcome comparisons by loading and displaying Medtronic(R) Framelink (TM) stereotactic plan with a PaCER electrode reconstructions (cf. Use-Case B).
 * **EXAMPLE_4.m** - Demonstrates PaCER operating in atlas space. Post OP CT and T1 linearly pre-registered to the template and electrodes plotted together with subcortical atlas segmentations. (cf. Use-Case D)
 * **EXAMPLE_5.m** - Demonstrates electrode reconstruction in **native** space with subsequent **transformation of the electrode objects to T1 space** (see example_1_2). The subcortical atlas from example 4 is now nonlinearly transformed and plotted with electrodes in T1 space. (cf. Use-Case D)
   * **EXAMPLE_5_1.m** - Continues EXAMPLE_5 by **integrating MRTrix fibertracking results** and plotting corticospinal tract streamlines along with electrodes and subcortical atlas segmentations.

### Questions
Feel free to drop a note to mail@andreashusch.de

### Acknowledgement
PaCER is packaged with some external toolboxes for convenience. Please see the "toolboxes" folder and the respective LICENSE files for details.


