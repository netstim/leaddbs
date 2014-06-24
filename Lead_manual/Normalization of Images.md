## Normalizing the Images

_Lead-DBS_ allows the user to automatically normalize the patient images into MNI space.

In this process, post-operative images are first coregistered with each other. These are then coregistered with the pre-operative images, in case they are available. Finally, all coregistrations are normalized into the MNI space.

_Lead-DBS_ has the following normalization protocols:
- _Schönecker 2009 linear threestep (Post-OP only):_
 This protocol...
- _Schönecker 2009 linear threestep (include Pre-OP):_
 This protocol...
- _Witt 2013 nonlinear:_
 This protocol...
- _Fuse CT with pre-operative MRI:_
 This protocol...

The user can select the protocol to run depending on the image files that are available for processing.


After loading the patient directory and choosing the imaging technique (see **Section 2**), check the option `[] Normalize` and press `Run`. During this process, the user will be asked to evaluate the coregistrations. Details on this step can be found in **Section 3.1**.


_Lead-DBS_ gives the user the option of also normalizing fiber tracking images into MNI space, when available. For processing of these images, the option `[] Normalize Fibers` must be checked.

The option `[] Check` lets the user look at the different views of the normalized images at the end of the process.

