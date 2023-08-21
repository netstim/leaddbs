LEAD-DBS
========

LEAD-DBS is ***NOT*** intended for clinical use!

## About Lead-DBS

LEAD-DBS is a MATLAB toolbox facilitating the:

- reconstruction of deep-brain-stimulation (DBS) electrodes in the human brain on basis of postoperative MRI and/or CT imaging
- the visualization of localization results in 2D/3D
- a group-analysis of DBS-electrode placement results and their effects on clinical results
- simulation of DBS stimulations (calculation of volume of activated tissue – VAT)
- diffusion tensor imaging (DTI) based connectivity estimates and fiber-tracking from the VAT to other brain regions (connectomic surgery)

## Installation

#### Prerequisites

- Recommended RAM size: 32GB or more
- MATLAB version: R2021a or later
- The following MATLAB toolboxes
  - MATLAB Image Processing Toolbox
  - MATLAB Signal Processing Toolbox
  - MATLAB Statistics and Machine Learning Toolbox
  - MATLAB Curve Fitting Toolbox (optional)
  - MATLAB Parallel Computing Toolbox (optional)
- The [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/download/) toolbox

#### Normal installation

Lead-DBS can be downloaded from our [website](https://www.lead-dbs.org) in fully functional form.

#### Development installation

Alternatively, especially in case you wish to modify and contribute to Lead-DBS, you can

- Make sure to meet the prerequisites
- Clone the Lead-DBS repository from [github](https://github.com/netstim/leaddbs.git).
- Download the necessary [data](https://www.lead-dbs.org/release/download.php?id=data_pcloud) and unzip it into the cloned git repository.

We’d love to implement your improvements into Lead-DBS – please contact us for direct push access to Github or feel free to add pull-requests to the Lead-DBS repository.

## Getting started

You can run Lead-DBS by typing "lead demo" into the Matlab prompt. This will open up the main GUI and a 3D viewer with an example patient.
But there's much more to explore. Head over to our [website](https://www.lead-dbs.org) to see a walkthrough tutorial, a manual, some more screenshots and other ressources. There's also a helpline in form of a Slack channel. We would love to hear from you.

## Questions

If you have questions/problems when using Lead-DBS, you can checkout our:

- Online [manual](https://netstim.gitbook.io/leaddbs/)
- Workthrough [videos](https://www.lead-dbs.org/helpsupport/knowledge-base/walkthrough-videos/)
- Knowledge [base](https://www.lead-dbs.org/helpsupport/knowledge-base/) (including [methods](https://www.lead-dbs.org/helpsupport/knowledge-base/lead-dbs-methods/), [cortical](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/cortical-atlas-parcellations-mni-space/)/[subcortical](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/atlases/) atlases, [connectomes](https://www.lead-dbs.org/helpsupport/knowledge-base/atlasesresources/normative-connectomes/), etc.)
- Support [forum](https://www.lead-dbs.org/?forum=lead-dbs-support-forum)
- [Auto-invite](https://www.lead-dbs.org/slack-invite) Slack [workspace](https://leadsuite.slack.com/)
