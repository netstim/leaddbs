## Writing 2D Images

After a successful reconstruction, you can obtain 2D images of each contact in its location within the brain tissue. To select this option, the box `[] Write out 2D` must be checked within the _Lead-DBS_ main window.

![Example of a 2D image](../images/2d_slices.png)
*Example of a 2D-slice export*

The 2D images are created using the localization settings obtained in the reconstruction and the manual correction. They help to better understand how the electrode contacts relate to the surrounding brain structures.

#### Features within the Images

Certain features can be changed to enhance the understanding of the structures surrounding the contacts:

- _Labels and color:_
 Specific brain areas are labeled and colored for easier identification. The labels can be turned off or on by checking the `Label` box. The color of the contour can also be selected from a pop-up window.
- _Bounding box:_
 Determines the final size of the 2D image.

Images are stored as `*.png` files within the patient folder, and are named according to the electrode and the plane they belong.
