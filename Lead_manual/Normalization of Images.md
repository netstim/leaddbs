## Normalizing the Images

_Lead-DBS_ allows the user to automatically normalize the patient images into MNI space.

In this process, post-operative images are first coregistered with each other. These are then coregistered with the pre-operative images, in case they are available. Finally, all coregistrations are normalized into the MNI space.

_Lead-DBS_ comes with the following built-in normalization protocols:
- _SPM DARTEL nonlinear (automatic):_
This protocol uses the Diffeomorphic Anatomical Registration Through Exponentiated Lie Algebra (DARTEL) approach supplied with SPM8 to normalize the preoperative MR-image directly to the ICBM template (in MNI space). The estimated DARTEL flowfields are then applied to the coregistered postoperative versions.
- _SPM DARTEL nonlinear:_
This protocol uses the Diffeomorphic Anatomical Registration Through Exponentiated Lie Algebra (DARTEL) approach supplied with SPM8 to normalize the preoperative MR-image directly to the ICBM template (in MNI space). The estimated DARTEL flowfields are then applied to the coregistered postoperative versions. In the process, you are asked by the program whether the image co-registration worked.
- _SPM Segment nonlinear (automatic):_
This protocol uses the SPM unified segmentation algorithm from SPM8 on the preoperative MR-image and applies the estimated deformation parameters to the coregistered postoperative versions.
- _SPM Segment nonlinear:_
This protocol uses the SPM unified segmentation algorithm from SPM8 on the preoperative MR-image and applies the estimated deformation parameters to the coregistered postoperative versions. In the process, you are asked by the program whether the image co-registration worked.
- _Schönecker 2009 linear threestep (Post-OP only):_
This protocol uses SPM to linearly coregister the postoperative images into MNI space in three consequtive steps, each focusing more on the subcortical target region. The last step spares the ventricles, which may largely vary in their subject-specific anatomy the most.
- _Schönecker 2009 linear threestep (include Pre-OP):_
This protocol uses SPM to linearly coregister the preoperative images into MNI space in three consequtive steps, each focusing more on the subcortical target region. The last step spares the ventricles, which may largely vary in their subject-specific anatomy the most. The estimated deformation matrix is then applied to coregistered postoperative images.
- _Fuse CT with pre-operative MRI:_
This protocol aims at fusing CT with MRI images and normalizing the images to MNI space based on parameters found for the MRI images. For now, this algorithm doesn't work robustly. To coregister MRI with CT images we propose the use of Slicer 3D or the Yale Bioimage Suite.

You can select the protocol to run depending on the image files that are available for processing.
You can easily include your own normalization procedure by writing a function "ea_normalize_yourprocedure(options)".

Add the following lines to the beginning of your code:

`if ischar(options) % return name of method.`
`    varargout{1}='My normalization technique';`
`    return;`
`end`

After loading the patient directory and choosing the imaging technique (see **Section 2**), check the option `[] Normalize` and press `Run`. During this process, the user will be asked to evaluate the coregistrations. Details on this step can be found in **Section 3.1**.

_Lead-DBS_ gives the user the option of also normalizing fiber tracking images into MNI space, when available. For processing of these images, the option `[] Normalize Fibers` must be checked.

The option `[] Check` lets the user look at the different views of the normalized images at the end of the process.

