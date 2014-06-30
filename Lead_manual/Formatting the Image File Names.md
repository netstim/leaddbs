## Renaming the Image Files

_Lead-DBS_ can process different image views. However, the image files must have a **specific name format** for them to be recognized. The naming scheme can be changed in the _**ea_prefs.m**_ file within the LEAD installation directory. If you use the built-in DICOM import function, you are asked to specify the type of each acquisition within a GUI and thus don't have to do the naming manually.

#### File Naming Format

**MR imaging** files within each folder must have the following format:
- Pre-operative images:
 - `pre_tra.nii` for transversal images
 - `pre_cor.nii` for coronal images
 - `pre_sag.nii` for sagittal images
 - `*.nii`for fiber tracking images

- Post-operative images:
 - `tra.nii` for transversal images
 - `cor.nii` for coronal images
 - `sag.nii` for sagittal images

**CT imaging** files within each folder must have the following format:

- Pre-operative MR images:
 - `pre.nii`
- Post-operative images:
 - `fusion.nii`

You can download a Nifti-Viewing software such as [MRIcron](http://www.mccauslandcenter.sc.edu/mricro/mricron/) to view the different image files and help in correct naming. [SPM8](http://www.fil.ion.ucl.ac.uk/spm/software/spm8/) needs to be installed (this program can also be used to visualize the images).

