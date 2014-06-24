## Checking the Coregistrations and Normalization

This process allows for better resulting images that will help to improve the localization of the electrodes.

Once the options are selected and the normalization process started, MATLAB will run the coregistration. During this process, the user will be shown coregistrations between different images and asked if they are precise. The user must evaluate the images and decide if the they are adequate.

**Adequate Coregistration:**
In case the coregistrations were adequately performed, the user should enter `y` into the Command Window. The present image will then be coregistered with the next available view.

This process happens up three times depending on the available image views. The order of coregistration follows:
- Transversal to coronal view
- Resulting image to sagittal view
- Resulting image to pre-operative transversal view
- Resulting image to the MNI space

The user is asked for evaluation between each item.

**Inadequate Coregistration:**
In case any of the coregistrations is not precise, the user should enter `n` into the Command Window. _Lead-DBS_ will then try again to coregister the images using a different routine. The user will be asked again to evaluate the images for adequate coregistration.
