*This project is in active development and does not have any regulatory approval*

# SlicerNetstim

This repository contains a collection of modules for 3D Slicer on the topic of deep brain stimulation and other applications&mdash;see the modules section below for a general overview.

## Modules provided in this extension

### Lead-OR
The Lead-OR module provides an interface for setting micro electrodes configuration in space. See the [Lead-OR page](./LeadOR/README.md) for installation instructions and how to load a sample dataset. See the [Lead-OR section in the Lead-DBS manual](https://netstim.gitbook.io/leaddbs/lead-or/imaging-setup) for the integration with electrophysiology recordings.

### WarpDrive
WarpDrive allows for manual interaction with the deformation fields that result from non linear registration. The idea is to be able to fix for miss-alignments between source and target images (and models). See the [WarpDrive page](./WarpDrive/README.md) for more information.

### Stereotactic Plan
This module outputs transforms that represent stereotactic planning trajectories from the planning settings. A routine to import plannings from Brainlab and ROSA is implemented. Output transforms are taken as an input to the Lead-OR module.

### Import Atlas
The import atlas module implements a routine to import [Lead-DBS atlases](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/atlases/). A [Lead-DBS installation](https://www.lead-dbs.org/download/) is needed so the path to the atlas can be provided. Atlases are imported into a model hierarchy in Slicer. [Atlases can also be loaded by drag-and-dropping the atlas_index.mat files](https://github.com/netstim/SlicerNetstim/pull/1).

### Import ACPC Autodetect
This module implements a routine to load AC,PC and MS coordinates from the `ACPC_autodetect.mat` file generated from Lead-DBS. The ACPC transform is also computed on load.

### FiducialRegistrationVariableRBF
This CLI module implements similar logic as [platimatch's landmark_warp](https://plastimatch.org/landmarks.html), but with variable radial basis functions. It used by WarpDrive to calculate transformations.

### CompositeToGridTransform
This CLI module creates a grid transform from a composite one. The logic is as the one in [Slicer's transforms module](https://github.com/Slicer/Slicer/blob/main/Modules/Loadable/Transforms/Logic/vtkSlicerTransformLogic.cxx#L561), but implemented in a CLI. This is useful to run it in the background and get progress indication, as used in WarpDrive.
## Illustrations

Such a scene can be represented when using Lead-DBS and modules of this extension. See the [Lead-OR page](./LeadOR/README.md) to load a sample dataset in Slicer.

![](LeadOR/Screenshots/Lead-OR_Scene.png?raw=true)

## References

- Oxenford, S., Roediger, J., Neudorfer, C., Milosevic, L., Güttler, C., Spindler, P., Vajkoczy, P., Neumann, W.-J., Kühn, A., & Horn, A. (2022). Lead-OR: A multimodal platform for deep brain stimulation surgery. *ELife*, 11, e72929. [https://doi.org/10.7554/eLife.72929](https://doi.org/10.7554/eLife.72929)

