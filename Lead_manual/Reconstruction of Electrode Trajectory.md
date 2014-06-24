## Reconstructing the Electrode Trajectory

Once normalized volumes are found within the chosen patient folder, the images are ready to be used for reconstruction of the electrode trajectories and manual correction, if needed.

#### Reconstructing Electrode Lead Trajectories
For the reconstruction step to take place, the user must check the box `[] Reconstruct`, choose the preferred parameters for obtaining the images, and press `Run`.

**IMPORTANT:**
Always pay attention to the checkboxes. Remember that _Lead-DBS_ runs all processes that are checked!

#### Reconstruction Parameters

To perform a reconstruction as precise as possible, _Lead-DBS_ uses the different planes to pinpoint the artifacts caused by the electrodes and calculates thereafter its trajectories through space. The user can also choose to reconstruct one or both hemispheres.

Several options are available to help in this process:

##### 1. Entry point for electrodes
The parameter `Entry point for electrodes` presents the user with following options:
```
- STN, GPi, or ViM
- Cg25
- Manual
```

An **automatic ** reconstruction will be performed if any of the first two options are chosen. The option `STN, GPi, or ViM` targets electrodes that have been implanted in patients with movement disorders. The option `Cg25`targets those in patients with depression.

The option `Manual` will require the user to pinpoint the entry points of the artifacts within the images. **Section 4.2** describes the details for this step.

##### 2. Axis

The parameter `Axis` determines the image planes that _Lead-DBS_ will use to locate the electrode. The following options are available:
```
- Use transversal image only
- Use transversal but smooth
- Use average of coronal and transversal, smoothed
```

##### 3. Mask window size
The default option is an **auto** window size. This has proven to give good results in the reconstruction.

However, numeric values (best results obtained from **5** to **15**) can be entered to fix the size of the mask. A smaller mask will avoid nearby structures that could interfere in the reconstruction step.



