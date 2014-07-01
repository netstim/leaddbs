# Reconstruction Statistics

LEAD-DBS exports statistics about your electrode reconstructions during
- 3D visualization
- Simulation of stimulations (calculation of Volume of Activated tissue and fiber connectivities).

During both processes, a file named `ea_stats.mat` is written out in MATLAB which can be used within your own analysis pipelines or together with the analysis tool `lead_group`. The structure of the variable file is as following:

- `conmat`: An _m_ × _n_ matrix showing the distance between the center of each contact (_m_ electrodes contacts numbered from K0 (right hemisphere) upwards) and the closest voxel within the _nth_ atlas. Corresponding atlas names can be found in the file, as well (see below). Please note that these values are never zero, even if the electrode contact resides within the atlas. This is due to the fact that both voxels and electrode contacts are represented as mathematical point and these are "never" exactly identical. Thus, to answer the question about whether a contact resides within a certain nucleus, you have to either define a distance threshold yourself, or try to work with one of the following two entries.
- `conmat_inside_vox`: A copy of `conmat`, where all entries are set to zero that show a smaller distance than the voxel resolution of the according atlas. This matrix can be used to see whether a contact resides surely within a certain atlas.
- `conmat_inside_hull`: Another copy of `conmat`, where all entries are set to zero for the condition that the center of an electrode contact resides within the concave hull of a certain nucleus. Please note that this is only valid for the _center_ of each contact. If the center of an electrode contact is close to the boarder of the nucleus, the contact might still reside "within" the nucleus without the entry of this matrix being zero. Taken together, for binary classifications, our advise is to analyze data with the pure `conmat` and defining your own reasonable threshold that tells you whether a contact is "within" or "outside" a certain nucleus.
- `patname`: The name/pseudonym of the analyzed patient.
- `atlases`: Information about the atlases used.
    - `names`: The names of the atlases. You can use this list for referencing all stats entries. As an example, entry (5,6) of `conmat` described above shows the distance between contact
