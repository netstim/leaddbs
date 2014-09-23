## Manual Review of Electrode Localization

When `[] Review` is checked and run, _Lead-DBS_ brings up a window where you can visualize and modify the position of the electrodes.

This step should **always** be done in order to obtain a meaningful reconstruction.

![Window for the Manual Correction](../images/manualcorrection.png)

**The goal of this step is to adjust the localization of the electrode and contact templates within the artifact, in order to place it as accurately as possible.**


In this window, the user is able to visualize the position of the electrode trajectory in the artifacts of the MR or CT images. _Lead-DBS_ places within these artifacts the built-in templates for the electrodes, represented as _straight lines_ that run along the artifacts. The contacts are represented as a _set of circles_ on the electrode lead lines (4 or 8 circles, depending on the electrode model chosen).

The **Manual Correction** window is divided into the following sections:
- Main view of the electrodes and contacts (center)
- Transversal views of electrode contacts (right column)
- Command buttons

In order to adjust the localization of the electrodes, two-dimensional planes are presented orthogonally to the lead trajectory. Changes are tracked within the main view and the transversal views. The main view presents the electrodes and contacts as a whole. The transversal views correspond to the position of the top and bottom contacts of each electrode. They are arranged in ascending order (top: `0`, second from top `3`, third from top `4`, and bottom `7`).

A set of commands are available in the manual correction window in order to adjust the settings of the electrode trajectories. Some also have keyboard shortcuts for easier access.

#### Changing Views

The main view displays the **different planes** that are selected by clicking on the command button or by using the following keyboard shortcuts:
- `A `: sets view from Anterior
- `P `: sets view from Posterior
- `R `: sets view from Right
- `L `: sets view from Left
- `X `: sets view from X direction
- `Y `: sets view from Y direction

**Brightness offset** and **contrast** are adjusted using their command buttons or the following keyboard shortcurts:

- `C`: decreases contrast
- `V`: increases contrast
- `B`: decreases offset
- `N`: increases offset

#### Adjusting the Electrode Height

If no contact is selected (scroll down to learn more on how to select a specific contact), the **arrow keys** will move the electrodes up or down (in the axial plane) as a whole within the artifact template in the following way:
- `Up` and `Down` arrow keys adjust the **right electrode**
- `Right` and `Left` arrow keys adjust the **left electrode**

Pressing the `Shift` key together with the **arrow** keys will move the electrodes in even greater steps.

Although the electrode templates already provide a specific **spacing** between each contact, these can be brought together or separated by using the following commands:

- `+`: increase spacing between contacts
- `-`: decrease spacing between contacts

#### Adjusting a Specific Contact's Position

For easier adjustment, specific contacts can be selected and their changes checked in their respective transversal views. Contacts are selected by clicking the command button or by pressing the keyboard shortcut:

- `0`: selects electrode 0 (right electrode)
- `3`: selects electrode 3 (right electrode)
- `4`: selects electrode 4 (left electrode)
- `7`: selects electrode 7 (left electrode)

When a specific contact is selected, the **arrow** keys will adjust the position within the _axial_ and _coronal_ planes. Changes can be tracked in the main view and the transversal views.

#### Saving the Changes

Once the manual correction is finished, the work can be saved by clicking the check on the command buttons or by using the keyboard shortcut:
- `space bar`: Finish manual corrections

Information will be saved as a Matlab file within the patient folder.
