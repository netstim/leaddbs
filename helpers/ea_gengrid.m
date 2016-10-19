function ea_gengrid(options)
if options.prefs.normalize.createwarpgrids
    if ~exist([options.root,options.patientname,filesep,'grid.nii'],'file')
        directory=[options.root,options.patientname,filesep];
        nii=ea_load_nii([directory,options.prefs.prenii_unnormalized]);
        nii.img(:)=0;
        spacing=10;
        nii.img(1:spacing:end,:,:)=255;
        nii.img(:,1:spacing:end,:)=255;
        nii.img(:,:,1:spacing:end)=255;
        nii.dt=[4 0];
        nii.fname=[directory,'grid.nii'];
        ea_write_nii(nii);
    end
end



