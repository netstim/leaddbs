## Manual Correction of Electrode Localization

When `[] Manual correction` is checked and run, _Lead-DBS_ brings up a window where the user can visualize and modify the position of the electrodes.
This should **always** be done to obtain a meaningful reconstruction.

![Image of the new window for the manual correction](http://www.andreas-horn.de/leaddbs/manualimages/manualcorrection.png)

**The goal of this step is to adjust the localization of the electrode and contact templates within the artifact, in order to place it as accurately as possible.**

#### 1. The Manual Correction Window

The window is divided into the following sections:
- Command buttons
- Main view of the electrodes and contacts (center)
- Transversal views of electrode contacts (right column)

##### Command buttons

A set of commands can be used in the manual correction window in order to adjust the settings of the electrode trajectories. Some also have keyboard shortcuts for easier access.

**Important:**
The **arrow keys** allow the user to adjust the position of the electrodes in the following way:
- Right and left arrow keys adjust the

When arrow keys are used, the whole electrode trajectory will be
If no contact is selected (scroll down to learn more on how to select a specific contact), the arrow keys will adjust the electrodes as a whole.


Electrodes can be adjusted as a whole or as specific contacts (0,3,4, or 7) by pressing t

Here is a list of the commands and their respective shortcuts (in most cases, you can increase the step size by simultaneously pressing the `shift` key):
- Decrease contrast `C`
- Increase contrast `V`
- Increase Offset `B`
- Decrease Offset `N`
- Select Electrode 0 `0`
- Select Electrode 3 `3`
- Select Electrode 4 `4`
- Select Electrode 7 `7`
- Increase spacing between contacts `+`
- Decrease spacing between contacts `-`
- Set view from Right `R`
- Set view from Left `L`
- Set view from Anterior `A`
- Set view from Posterior `P`
- Set view from X direction `X`
- Set view from Y direction `Y`
- Finish manual corrections `space bar`


##### Main view

In this window, the user is able to visualize the position of the electrode trajectory in the artifacts of the MR or CT images. _Lead-DBS_ places within these artifacts the built-in templates for the electrodes, represented as _straight lines_ that run along the artifacts. The contacts are represented as a _set of circles_ on the electrode lead lines (4 or 8 circles, depending on the electrode model chosen).

In order to adjust the height of the electrodes, two-dimensional planes are presented orthogonally to the trajectory.

Using the command buttons (or their respective keys), you can change the view of the scene, change the brightness and contrast settings, adjust the position of the electrodes,
