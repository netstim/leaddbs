## Checking the Coregistrations and Normalization

This process allows for better resulting images that will help to improve the localization of the electrodes.

Once the options are selected and the normalization process started, MATLAB will run the coregistration. Depending on the chosen normalization process, you will be shown coregistrations between different images and asked if they are precise. You should then evaluate the images and decide if the they are adequately coregistered.

**Adequate Coregistration:**
In case the coregistrations were adequately performed, you should enter `y` into the Command Window. The present image will then be coregistered with the next available view.

This process happens up to three times depending on the available image views. The order of coregistration follows:
- Transversal to coronal view
- Resulting image to sagittal view
- Resulting image to pre-operative transversal view
- Resulting image to the MNI space

You are asked for evaluation between each item.

**Inadequate Coregistration:**
In case any of the coregistrations is not precise, you should enter `n` into the Command Window. _Lead-DBS_ will try again to coregister the images using a different optimization metric and you will be asked again to evaluate the images. The toolbox attempts to run the coregistrations up to **four** times (that is, entering `n` 3 times).
