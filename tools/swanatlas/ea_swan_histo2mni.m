function ea_swan_histo2mni(swanpath)

indir=uigetdir(fullfile(swanpath,'histology','masks'));

outfile=ea_swan_imgs2nii(indir, [0.3,0.3,0.3]);

options=ea_getptopts(fullfile(swanpath,'derivatives','leaddbs','sub-SWANatlas'));

[pth,fn,ext]=fileparts(outfile);
ea_apply_normalization_tofile(options, outfile, fullfile(pth,[fn,'_template',ext]), 0, 1, fullfile(ea_space,'t1.nii'));

disp(['Normalized segmentation exported to: ',fullfile(pth,[fn,'_template',ext]),'.']);

