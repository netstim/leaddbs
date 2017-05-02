% generate atlas

aID = fopen('inia19-NeuroMaps_selection.txt');
atlas_lgnd=textscan(aID,'%d %s %d %d %d %d');
mkdir('../atlases/INIA19 (Rohlfing 2012)');
mkdir('../atlases/INIA19 (Rohlfing 2012)/lh');
mkdir('../atlases/INIA19 (Rohlfing 2012)/rh');
atl=ea_load_nii('inia19-NeuroMaps.nii');
atl.img=round(atl.img);

% symmetrize hemispheres
todelete=[];
for reg=1:length(atlas_lgnd{1})
   if ~ismember(atlas_lgnd{1}(reg)+1000,atlas_lgnd{1}) && ~ismember(atlas_lgnd{1}(reg)-1000,atlas_lgnd{1})
       todelete=[todelete,reg];
   end
end
atlas_lgnd{1}(todelete)=[];
atlas_lgnd{2}(todelete)=[];

for reg=1:length(atlas_lgnd{1})

    snii=atl;
    snii.img(snii.img~=atlas_lgnd{1}(reg))=0;
    snii.img=logical(snii.img);

    if strcmp(atlas_lgnd{2}{reg}(1),'r')
        odir='../atlases/INIA19 (Rohlfing 2012)/rh';
    elseif strcmp(atlas_lgnd{2}{reg}(1),'l')
        odir='../atlases/INIA19 (Rohlfing 2012)/lh';
    else
        ea_error('Something is not right');
    end

    if any(snii.img(:))

        snii.fname=[odir,filesep,atlas_lgnd{2}{reg}(3:end),'.nii'];
        spm_write_vol(snii,snii.img);
        [pth,fn,ext]=fileparts(snii.fname);

        % warp to MNI:
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.matname = {'inia19_to_mni_hires_sn.mat'};
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.vox = [0.5 0.5 0.5];
        matlabbatch{1}.spm.util.defs.comp{1}.sn2def.bb = [NaN NaN NaN
            NaN NaN NaN];
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {snii.fname};
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {pth};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = '';

        spm_jobman('run',{matlabbatch});
        clear matlabbatch

        ea_crop_nii([pth,filesep,'w',fn,ext]);
        movefile([pth,filesep,'w',fn,ext],[pth,filesep,fn,ext]);
    end
end
