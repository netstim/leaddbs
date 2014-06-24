## Lead-DBS Main Window and Loading Images
The main window of Lead-DBS is divided into 5 sections:
- Selection of patient directory and imaging technique
- Normalization
- Reconstruction
- Manual correction
- Visualization

The button `Run` in the bottom right runs **ALL** processes that are checked within the main window.

**Notice**: If the user wants to run a specific process(es), he must tend that only the desired checkboxes are selected!

![Lead-DBS Main Window]()

#### 1. Selecting the file directory

By clicking on `Choose Patient Directory` on the top, the user can browse and select the location of the image folder.

For _Lead-DBS_ to analyze the images, it needs at least the following views:
- for MR images, a post-operative transversal view
- for CT images, a post-operative transversal view

The user may also have more views available, in which case the process of electrode localization can be made easier.

**IMPORTANT:**
Images must be in the `*.nii` file format and must have a **_specific naming format_**. For details on how to format the file names, please refer to **Section 2.1**.

#### 2. Selecting the imaging technique

After the image directory has been chosen, the user must select the imaging technology in the dropdown menu on the top right of the window. `MR` imaging is the default setting.
