function ea_coreg_all_mri(options,usebrainmask)

% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

modalities = fieldnames(options.subj.preopAnat);
presentfiles = cellfun(@(modality) options.subj.preopAnat.(modality).preproc, modalities, 'UniformOutput', false);

for coregfi=2:length(presentfiles)
    if ~ea_coreglocked(options,presentfiles{coregfi}) % file has already been locked and approved by used
        out_file = options.subj.coreg.anat.preop.(modalities{coregfi});
        ea_coreg2images(options, presentfiles{coregfi}, presentfiles{1}, out_file);

        % reslice images if needed
        V1 = ea_open_vol(presentfiles{1});
        V2 = ea_open_vol(out_file);
        if ~isequal(V1.mat,V2.mat)
            ea_conformspaceto(presentfiles{1}, out_file, 1);
        end
        % better slab support:
        nii=ea_load_nii(V2.fname);
        nii.img(abs(nii.img)<0.0001)=0;
        ea_write_nii(nii);
    end
end
