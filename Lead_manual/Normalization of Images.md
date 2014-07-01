## Normalizing the Images

_Lead-DBS_ allows you to automatically normalize the patient images into MNI space. For all normalizations, the [_ICBM 2009b Nonlinear Asymmetric_](http://www.bic.mni.mcgill.ca/ServicesAtlases/ICBM152NLin2009) is used. This template can be found in your installation under
`lead_dbs/templates/mnihires.nii`.

There are several approaches to achieve a precise normalization of pre- and post-operative images. In most of the routines that come pre-installed with LEAD, post-operative images are first coregistered with each other and then coregistered to the pre-operative images. Finally, all coregistrations are normalized into the MNI space based on the transformation parameters found when co-registering the pre-operative image to the template.

_Lead-DBS_ comes with the following built-in normalization protocols:

- _SPM DARTEL nonlinear:_

  This protocol uses the Diffeomorphic Anatomical Registration Through Exponentiated Lie Algebra (DARTEL) approach supplied with SPM8 to normalize the preoperative MR-image directly to the ICBM template (in MNI space). The estimated DARTEL flowfields are then applied to the coregistered postoperative versions.

  - _SPM New Segment nonlinear:_

  This protocol uses the SPM "New Segment" approach to segment and normalize the pre-operative image to the ICBM template (in MNI space). The estimated deformation fields are then applied to the coregistered postoperative versions. LEAD uses a slightly modified version of the New Segment approach in that it uses a higher spatial resolution of the warps. This leads to a higher processing time.

- _SPM Segment nonlinear:_

  This protocol uses the SPM unified segmentation algorithm from SPM8 on the preoperative MR-image and applies the estimated deformation parameters to the coregistered postoperative versions.

- _Schönecker 2009 linear threestep (Post-OP only):_

  This protocol uses SPM to linearly coregister the postoperative images into MNI space in three consecutive steps, each focusing more on the subcortical target region. The last step spares the ventricles, which may largely vary in the subject-specific anatomy. This is the only normalization routine that can handle the situation where you don't have pre-operative images and still gives very precise results (Applying _non-linear_ deformations to the post-operative images directly is not a good choice in our experience, since they lead to curly dbs-trajectories and thus render the results useless.

- _Schönecker 2009 linear threestep (include Pre-OP):_

  This protocol uses SPM to linearly coregister the preoperative images into MNI space in three consecutive steps, each focusing more on the subcortical target region. The last step spares the ventricles, which may largely vary in the subject-specific anatomy. The estimated deformation matrix is then applied to the coregistered postoperative images.

- _Fuse CT with pre-operative MRI:_

  This protocol aims at fusing CT with MRI images and normalizing the images to MNI space, based on parameters found for the MRI images. For now, this algorithm doesn't work robustly. To coregister MRI with CT images, we propose the use of [Slicer 3D](http://www.slicer.org/) or the [Yale BioImage Suite](http://medicine.yale.edu/bioimaging/suite/).

You can select the protocol to run depending on the image files that are available for processing.

You can easily include your own normalization procedure by writing a function "_ea_normalize_yourprocedure(options)_". Simply add the following lines to the beginning of your code:

```
if ischar(options) % return name of method.
    varargout{1}='My normalization technique';
    return;
end
```

After loading the patient directory and choosing the imaging technique (see **Section 2**), check the option `[] Normalize` and press `Run`. During the process, you will be asked to evaluate the coregistrations. Details on this step can be found in **Section 3.1**.

When available, _Lead-DBS_ also gives you the option of normalizing fiber tracking images into MNI space. For processing of these images, the option `[] Normalize Fibers` must be checked.

The option `[] Check` lets the user look at the different views of the normalized images at the end of the process.

