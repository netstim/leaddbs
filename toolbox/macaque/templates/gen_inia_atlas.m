% generate atlas

aID = fopen('inia19-NeuroMaps.txt');
atlas_lgnd=textscan(aID,'%d %s %d %d %d %d');
mkdir('../atlases/inia19_full');
mkdir('../atlases/inia19_full/lh');
mkdir('../atlases/inia19_full/rh');
atl=ea_load_nii('inia19-NeuroMaps.nii');
atl.img=round(atl.img);
for reg=1:length(atlas_lgnd{1})
    
    snii=atl;
    snii.img(snii.img~=reg)=0;
    
    if strcmp(atlas_lgnd{2}{reg}(1),'r')
        odir='../atlases/inia19_full/rh';
    elseif strcmp(atlas_lgnd{2}{reg}(1),'l')
        odir='../atlases/inia19_full/lh';
    else
        ea_error('Something is not right');
    end

    if any(snii.img(:))
        
        snii.fname=[odir,filesep,atlas_lgnd{2}{reg}(3:end),'.nii'];
        spm_write_vol(snii,snii.img);
        
        % warp to MNI:
        [pth,fn,ext]=fileparts(snii.fname);
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
        
        cfg_util('run',{matlabbatch});
        clear matlabbatch
        
        ea_crop_nii([pth,filesep,'w',fn,ext]);
        movefile([pth,filesep,'w',fn,ext],[pth,filesep,fn,ext]);
    end
end
