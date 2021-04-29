%write your .mat into nifti
ref = [ea_space,'t1.nii'];
refnii = ea_load_nii(ref);
refnii.img = '/your/.matfile';
refnii.fname = '/full/path/of/where/you/want/to/save/nii/file';
ea_write_nii(refnii);
result_fig = ea_mnifigure;
curr_ax = result_fig.CurrentAxes;
nifti_roi = ea_roi('path/.nii/file');
patch(curr_ax,'Faces',nifti_roi.fv.faces,'Vertices',nifti_roi.fv.vertices,'FaceColor','r','EdgeColor','r');
