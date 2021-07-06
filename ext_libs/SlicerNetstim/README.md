# SlicerNetstim

This repository contains a collection of modules for 3D Slicer on the topic of deep brain stimulation and other applications. See the modules section below for a general overview.

Part of this work was part of the [NA-MIC Project Week](https://projectweek.na-mic.org/PW35_2021_Virtual/) and is still under development.

## Modules

### AlphaOmega

This module implements an interface with the NeuroOmega device from the AlphaOmega company. It can query information about micro electrodes in real time. This module is disabled in the Slicer Extension given that uses a proprietary SDK. Work is being done to include it in a custom Slicer build.

### Import Atlas

The import atlas module implements a routine to import [Lead-DBS atlases](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/atlases/). A [Lead-DBS installation](https://www.lead-dbs.org/download/) is needed so the path to the atlas can be provided. Atlases are imported into a model hierarchy in Slicer. [Atlases can also be loaded by drag-and-dropping the atlas_index.mat files](https://github.com/netstim/SlicerNetstim/pull/1).

### LeadOR

The LeadOR module provides an interface for setting micro electrodes configuration in space. It takes trajectories created from the AlphaOmega module and the output transform from the StereotacticPlan to create a live update of the scene.

### StereotacticPlan

This module takes as input stereotactic frame coordinates (currently using Leksell) and creates a transform representing the trajectory in Slicer. An PDF import routine is implemented from Brianlab planning files.

### WarpDrive

The WarpDrive module is not necessarily DBS related. WarpDrive allows for manual interaction with the deformation fields that result from non linear registration. The idea is to be able to fix for miss-alignments between source and target images (and models). See the [WarpDrive page](./WarpDrive/README.md) for more information.

## Illustrations

Such a scene can be represented when using Lead-DBS and modules of this extension.

![](Documentation/Scene.png?raw=true)
![](Documentation/Screenshot.png?raw=true)