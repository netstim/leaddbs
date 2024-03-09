function ea_usp_histo2mni(uspath)

indir=uigetdir(fullfile(uspath,'histology','masks'));

outfile=ea_usp_imgs2nii(indir, [0.3,0.3,0.3]);

options=ea_getptopts(fullfile(uspath,'derivatives','leaddbs','sub-USPatlas'));

[pth,fn,ext]=fileparts(outfile);
ea_apply_normalization_tofile(options, outfile, fullfile(pth,[fn,'_template',ext]), 0, 1, fullfile(ea_space,'t1.nii'));

disp(['Normalized segmentation exported to: ',fullfile(pth,[fn,'_template',ext]),'.']);

