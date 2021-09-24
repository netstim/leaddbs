[back to main page](https://github.com/netstim/SlicerNetstim)

# Quick Start Guide

## Slicer Installation 

To install the SlicerNetstim Extension it is first necessary to install Slicer. Slicer can be downloaded from their [downloads page](https://download.slicer.org/). We recommend to use the 4.11 stable version, 4.13 nightly (and upcoming) versions should work as well.

Once Slicer is installed, open it and go to the Extensions Manager clicking it's icon in the toolbar (see picture below).

![](ExtensionManager.png?raw=true)

## SlicerNetstim Installation

Install the SlicerNetstim Extension from the Extension Manager. It is possible to use the search function to filter the extensions by name or look for it within the available extensions.

![](SlicerNetstimInstall.png?raw=true)

Once installed, restart Slicer.

Great! Now, the Netstim category should show under the Module Selector.

![](SlicerNetstimInstalled.png?raw=true)

## Download Sample Data

Before diving into Lead-OR let's download a sample dataset. For this, first switch to the Sample Data module.

![](SampleData.png?raw=true)

Scroll down and select the STN Planning dataset under the Lead-OR section. This will download the dataset and can take a minute or two.

![](SampleDataDownload.png?raw=true)

This dataset contains T1, T2 and fgatir sequences together with the DISTAL atlas warped into patient space and a Planning transformation. You can see the data in the Slicer Scene by switching to the Data module.

![](Data.png?raw=true)


## Lead-OR

Now back to the module selector and switch to the Lead-OR module. Here, set the Planning transform with the available planning and the Distance to Target with the respective transform.

![](Lead-OR.png?raw=true)

By selecting the micro electrode controllers you can see the micro electrodes appear in the 3D view and the slice displays as well.

By sliding the Distance to Target slider, the micro electrodes translate in space.

![](MicroElectrode.png?raw=true)
