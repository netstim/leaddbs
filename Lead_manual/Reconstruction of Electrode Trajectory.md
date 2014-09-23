## Reconstructing the Electrode Trajectory

Once normalized volumes are found within the chosen patient folder, the images are ready to be used for reconstruction of the electrode trajectories and for manual correction.

#### Reconstructing Electrode Lead Trajectories
For the reconstruction step to take place, check the box `[] Reconstruct`, choose the preferred parameters for obtaining the images, and press `Run`.

By default, _Lead-DBS_ performs an automatic search of the artifacts caused by the electrode leads. This can be modified using the following parameters.
- _Entry point:_ set `STN, GPi, or ViM`, `Cg25` or `Manual`. If set to _Manual_, you are prompted to select the starting point of the artifact for each side.
- _Axis:_ Several “preprocessing”-options can be selected to optimize the imaging data you want to use for the reconstruction here.
- _Mask window size:_ either set a concrete value (within the range of 5 to 20 might improve results) or set the value to “auto”. This option is only important for MR images, in CT imaging, a large mask size (e.g. 20) can be used but “auto” should work as well. If performing reconstructions based on MR images, please use a larger mask value if there is a lot of edema around the electrode and a smaller value if the images are rather noisy. Smaller values will prevent the algorithm to “lose track” but might stop to soon if the artifacts are rather large.

In most cases, an automatic reconstruction will render an adequate result. If this is not the case, you can improve the process by choosing different options within the parameters `Axis`or `Mask window size` (see above and below for more detail).


#### Reconstruction Parameters in more detail

To perform a reconstruction as precise as possible, _Lead-DBS_ uses different MR-acquisitions (if available) to pinpoint the artifacts caused by the electrodes. You can choose to reconstruct one or both hemispheres (`[] LH`and `[] RH` checkboxes).

Several options are available to help in this process:

##### 1. Entry point for electrodes
The parameter `Entry point for electrodes` presents following options:
```
- STN, GPi or ViM
- Cg25
- Manual
```

An **automatic ** reconstruction will be performed if any of the first two options are chosen. The option `STN, GPi, or ViM` targets electrodes that have been implanted in patients with movement disorders. The option `Cg25` targets those in patients with depression.

The option `Manual` will require you to pinpoint the entry points of the artifacts within the image slices. This option should in theory work with all DBS electrodes (wherever they might be implanted). Please select the right electrode when prompted the first time and the left one when prompted the second time. The red square shows the area, LEAD DBS would scan in the `STN, GPi or ViM` mode and helps you to identify the right hemisphere in each step.

##### 2. Axis

The parameter `Axis` determines the images and possible small preprocessing steps that _Lead-DBS_ will use to locate the electrode. The following options should be preferred:
```
- Use transversal image only
- Use transversal but smooth
- Use average of coronal and transversal, smoothed
```

##### 3. Mask window size
The default option is an _automatic_ window size (enter “auto” into the field). This has proven to often give good results in the reconstruction.

However, numeric values (best results obtained with values between **5** to **15** or **20**) can be entered to fix the size of the mask. A smaller mask will avoid nearby structures that could interfere in the reconstruction step. If the image shows large artifacts, e.g. due to local edema, a larger mask should be chosen (e.g. enter `15` instead of `auto`). If the image is noisy, a smaller mask (e.g. `5` or `7`) might be of better use.



