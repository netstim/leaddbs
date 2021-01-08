function ea_normsubcorticalsegm(options)

directory=[options.root,options.patientname,filesep];

if exist([directory,'mmnormalized.mat'],'file')
    return
end

mkdir([directory,'atlases',filesep,'mni',filesep]);
mniatldir=[directory,'atlases',filesep,'mni',filesep,options.patientname,filesep];
natatldir=[directory,'atlases',filesep,'native',filesep,options.patientname,filesep];
mkdir(mniatldir);
mkdir([mniatldir,'lh']);
mkdir([mniatldir,'rh']);
whichnormmethod=ea_whichnormmethod(directory);

srcs={'Pallidum','Ruber','STN'};
for iside=1:length(options.sides)
    switch options.sides(iside)
        case 1
            sidec='rh';
        case 2
            sidec='lh';
    end
    for src=1:length(srcs)

        switch whichnormmethod
            case ea_getantsnormfuns
                ea_ants_apply_transforms(options,{[natatldir,sidec,filesep,srcs{src},'.nii']},{[mniatldir,sidec,filesep,srcs{src},'.nii']},0,[ea_space,'bb.nii']);
            case ea_getfslnormfuns
                ea_fsl_apply_normalization(options,{[natatldir,sidec,filesep,srcs{src},'.nii']},{[mniatldir,sidec,filesep,srcs{src},'.nii']},0,[ea_space,'bb.nii']);
            otherwise
                matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.fnames = {[natatldir,sidec,filesep,srcs{src},'.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
                matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {[natatldir,sidec,filesep]};
                matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[ea_space,'bb.nii']};
                matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
                matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = [0 0 0];
                spm_jobman('run',{matlabbatch});
                clear matlabbatch
                movefile([mniatldir,sidec,filesep,'w',srcs{src},'.nii'],[mniatldir,sidec,filesep,srcs{src},'.nii']);
        end
        ea_crop_nii([mniatldir,sidec,filesep,srcs{src},'.nii']);

    end
end

mmnormalized=1;
save([directory,'mmnormalized.mat'],'mmnormalized');


