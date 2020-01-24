LEAD-DBS
========

LEAD-DBS is *NOT* intended for clinical use!

### About Lead-DBS

LEAD-DBS is a MATLAB-toolbox facilitating the: 

- reconstruction of deep-brain-stimulation (DBS) electrodes in the human brain on basis of postoperative MRI and/or CT imaging
- the visualization of localization results in 2D/3D
- a group-analysis of DBS-electrode placement results and their effects on clinical results
- simulation of DBS stimulations (calculation of volume of activated tissue – VAT)
- diffusion tensor imaging (DTI) based connectivity estimates and fiber-tracking from the VAT to other brain regions (connectomic surgery)

LEAD-DBS builds on SPM8/12, especially regarding warping and segmentation procedures.

### Installation

Usually, Lead-DBS can be downloaded from our website (www.lead-dbs.org) in fully functional form.
Alternatively, especially in case you wish to modify and contribute to Lead-DBS, you can

- Clone the Lead-DBS repository from here.
- Download the [necessary datafiles](http://www.lead-dbs.org/release/download.php?id=data) using this link and unzip the downloaded folder into the cloned git repository.
- We’d love to implement your improvements into Lead-DBS – please contact us for direct push access to Github or feel free to add pull-requests to the Lead-DBS repository.

### Getting started

You can run Lead-DBS by typing "lead demo" into the Matlab prompt. This will open up the main GUI and a 3D viewer with an example patient.
But there's much more to explore. Head over to https://www.lead-dbs.org/ to see a walkthrough tutorial, a manual, some more screenshots and other ressources. There's also a helpline in form of a Slack channel. We would love to hear from you.
