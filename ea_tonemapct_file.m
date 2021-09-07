function ea_tonemapct_file(options,nativemni)
ea_dispt('Tonemapping CT...');
if ~exist('nativemni','var')
    nativemni='native';
else
    if isempty(nativemni)
            nativemni='native';
    end
end

switch nativemni
    case 'native'
        if isfile(options.subj.postopAnat.CT.coreg)
            ct = ea_load_nii(options.subj.postopAnat.CT.coreg);
            ct.fname = options.subj.postopAnat.CT.coregTonemap;
            ct.img = ea_tonemap_ct(ct.img);
            ct.dt = [16,0];
            ea_write_nii(ct);
        else
            fprintf('Coregistered CT image not present. Skipping tonemapping.\n');
        end
    case 'mni'
        if isfile(options.subj.postopAnat.CT.norm)
            ct = ea_load_nii(options.subj.postopAnat.CT.norm);
            ct.fname = options.subj.postopAnat.CT.normTonemap;
            ct.img = ea_tonemap_ct(ct.img);
            ea_write_nii(ct);
        else
            fprintf('Normalized CT image not present. Skipping tonemapping.\n');
        end
end

ea_dispt('');
