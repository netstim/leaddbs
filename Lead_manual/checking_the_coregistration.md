## Checking the Quality of Normalizations and Co-registrations

This process allows for better resulting images that will help to improve the localization of the electrodes.

### Normalization

It is very crucial to manually check the normalization quality of your images before starting the reconstruction process. No matter how you normalized your images, you can always set the `Check` checkbox below the `Normalization` checkbox to display slice images of your normalized images that are overlayed with a wire-frame like structure that has been generated from the ICBM template. Only if these results look very precise, i.e. you can see that the wire-frames from the template match the anatomy of the patients normalized MR images, you should advance to the reconstruction step.

### Coregistrations between pre- and postop MR images

Once the options are selected and the normalization process started, MATLAB will run the coregistration. Depending on the chosen normalization process, you will be shown coregistrations between different images and asked if they are precise, if you set `prefs.normalize.coreg` in the `ea_prefs.m` file to _'manual'_. You can then evaluate the images and decide if the they are adequately coregistered. In most cases, you can skip this step completely, i.e. leave the setting in `ea_prefs.m` set to _'auto'_.

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
