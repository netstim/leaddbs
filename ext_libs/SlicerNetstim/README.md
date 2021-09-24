*This project is in active development and does not have any regulatory approval*

# SlicerNetstim

This repository contains a collection of modules for 3D Slicer on the topic of deep brain stimulation and other applications. See the modules section below for a general overview and the [quick start guide](./Documentation/QuickStartGuide/README.md) to install and load sample data in Slicer.

## Quick Start

Follow the [quick start guide](./Documentation/QuickStartGuide/README.md) to install the software and load a sample dataset with STN DBS planning in Slicer.

## Modules

- AlphaOmega

This module implements an interface with the NeuroOmega device from the AlphaOmega company. It can query information about micro electrodes in real time. This module is disabled in the Slicer Extension given that uses a proprietary SDK. To build this extension follow the [Slicer development guidelines](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html).

- Import Atlas

The import atlas module implements a routine to import [Lead-DBS atlases](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/atlases/). A [Lead-DBS installation](https://www.lead-dbs.org/download/) is needed so the path to the atlas can be provided. Atlases are imported into a model hierarchy in Slicer. [Atlases can also be loaded by drag-and-dropping the atlas_index.mat files](https://github.com/netstim/SlicerNetstim/pull/1).

- Lead-OR

The Lead-OR module provides an interface for setting micro electrodes configuration in space. It takes trajectories created from the AlphaOmega module and the output transform from the StereotacticPlan to create a live update of the scene.

- StereotacticPlan

This module takes as input stereotactic frame coordinates (currently using Leksell) and creates a transform representing the trajectory in Slicer. An PDF import routine is implemented from Brianlab planning files.

- WarpDrive

WarpDrive allows for manual interaction with the deformation fields that result from non linear registration. The idea is to be able to fix for miss-alignments between source and target images (and models). See the [WarpDrive page](./WarpDrive/README.md) for more information.

## Illustrations

Such a scene can be represented when using Lead-DBS and modules of this extension. See the [quick start guide](./Documentation/QuickStartGuide/README.md) to load a sample dataset in Slicer.

![](Documentation/Scene.png?raw=true)
![](Documentation/Screenshot.png?raw=true)