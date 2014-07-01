## Reconstruction Statistics

LEAD-DBS exports statistics about your electrode reconstructions during
- 3D visualization
- Simulation of stimulations (calculation of Volume of Activated tissue and fiber connectivities).

During both processes, a file named `ea_stats.mat` is written out in MATLAB which can be used within your own analysis pipelines or together with the group analysis GUI `lead_group`. The structure of the variable file is as following:

- `conmat`: An _m_ × _n_ matrix showing the distance between the center of each contact (_m_ electrodes contacts numbered from K0 (right hemisphere) upwards) and the closest voxel within the _nth_ atlas. Corresponding atlas names can be found in the file, as well (see below). Please note that these values are never zero, even if the electrode contact resides within the atlas. This is due to the fact that both voxels and electrode contacts are represented as mathematical point and these are "never" exactly identical. Thus, to answer the question about whether a contact resides within a certain nucleus, you have to either define a distance threshold yourself, or try to work with one of the following two entries.
- `conmat_inside_vox`: A copy of `conmat`, where all entries are set to zero that show a smaller distance than the voxel resolution of the according atlas. This matrix can be used to see whether a contact resides surely within a certain atlas.
- `conmat_inside_hull`: Another copy of `conmat`, where all entries are set to zero for the condition that the center of an electrode contact resides within the concave hull of a certain nucleus. Please note that this is only valid for the _center_ of each contact. If the center of an electrode contact is close to the boarder of the nucleus, the contact might still reside "within" the nucleus without the entry of this matrix being zero. Taken together, for binary classifications, our advise is to analyze data with the pure `conmat` and defining your own reasonable threshold that tells you whether a contact is "within" or "outside" a certain nucleus.
- `patname`: The name/pseudonym of the analyzed patient.
- `atlases`: Information about the atlases used.
    - `names`: The names of the atlases. You can use this list for referencing all stats entries. As an example, entry (5,6) of `conmat` described above shows the distance between contact
    - `types`: The types of atlas. Types can be the following:
        - _1_: Right hemispheric atlas, i.e. one that is _only_ present in the _atlases/youratlasset/rh_ directory.
        - _2_: Left hemispheric atlas, i.e. one that is _only_ present in the _atlases/youratlasset/lh_ directory.
        - _3_: An atlas that is present for both hemispheres and is composed of two files with the exact same name, one lying in the _atlases/youratlasset/lh_ and one in the _atlases/youratlasset/rh_ directory.
        - _4_: An atlas that is present for both hemispheres but is read-in only from one file residing within the _atlases/youratlasset/mixed_ directory. When processed, these atlases are divided along the _x = 0 mm_ (midsaggital) plane into two compartments.
        - _5_: An atlas that shows a midline structure, i.e. one that is not being divided along the _x = 0 mm_ plane. These atlases can be used for structures around the midsaggital plane, such as the subcallosal cortex or Broca's area 25.
        - `colormap`: A colormap used for visualization. You can change the colormap from within the `lead` GUI by pressing on the colormap button, selecting a colormap and re-running the 3D visualization step.
    - `electrodes`: A structure with information about the reconstruction. This is the same information as can be found in the `ea_reconstruction.mat` file within the patient folder, namely:
        - `coords_mm`: A (usually) 1 × 2 cell filled with n × 3 matrices denoting the coordinates of each electrode contact in millimeters (within MNI space). The first cell accounts for the right electrode, the second one for the left electrode. Within each matrix, the first row accounts for the dorsalmost contact. For example, _ea_stats.electrodes.coords_mm{2}(3,:)_ accounts for the coordinates of the third-dorsalmost contact on the left hemisphere, often referred to as `K10` (or sometimes as `K6`).
        - `trajectory`: A (usually) 1 × 2 cell filled with n × 3 matrices, each describing a coordinate list of the reconstructed trajectory. Points on these line should be equidistant and on a straight line.
        - `name`: The name/pseudonym of the patient.
    - `vat`: Information about volumes of activated tissue (VAT). Each struct in here has the following fields:
        - `U`: The Voltage applied to build this VAT in V.
        - `Im`: The Impedance measured while stimulating, in Ω.
        - `Contact`: The electrode contact where this VAT found its origin. An entry [1 2] would refer to right (1) contact number 2, i.e. the second-dorsalmost right contact, usually referred to as `K1`. [2 2] would be its left counterpart `K9`/`K5`.
        - `Side`: Whether the contact was left or right (redundant with the first entry in the `Contact` field).
        - `AtlasIntersection`: The 1 × n matrix in this field estimates, how many cubic millimeters of tissue of each nucleus got stimulated by the generated VAT (n := number of atlases in the set). If the second entry in this matrix amounts to a value of 8.7, this means that according to the simulation, 8.7 cubic millimeters of the second atlas (as labeled by _ea_stats.atlases.names{2}) was stimulated by this particular VAT.
    - `ft`: A structure accounting for fibertracking-related statistics.
        - `fibercounts`: A n × 1 matrix accounting for how many fibers connected the VAT to the n structures within the particular labeling atlas.
        - `labels`: A n × 1 cell with the names of the labeling atlas.
        - `vatsused`: refers to the VATs used for this fiber-tracking analysis. If you have more than one struct of `ft`, this can be used to keep track on which analysis was which.
    - `vatanalyses`: A meta-structure keeping track on the analyses performed. For example, if this structure has three entries, three analyses have been performed. In theory, this could also be used to keep track of stimulation settings that may change over time for a certain patient. In `lead_group`, the _Analysis_ listbox refers to the last entry in `vatanalyses`, if it is set to 0. It refers to the second-last entry if it is set to -1, and so forth. This way, more than one VAT/Fibertracking analyses can be stored for each patient and analysed seperately with `lead_group`.
        - `vatsused`: Which entries of the _ea_stats.vat_ structure have been used in this particular analysis.
        - `fibersused`: Which entries of the _ea_stats.ft_ structure have been used in this particular analysis.
