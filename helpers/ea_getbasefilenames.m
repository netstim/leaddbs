function leadfiles = ea_getbasefilenames(options)
% Return names of the files pre-defined in lead

if exist('options','var')
    if ~isempty(options)
        prefs=ea_prefs(options.patientname);
    else
        prefs=ea_prefs('');
    end
else
    prefs=ea_prefs('');
end


corefiles={
    prefs.prenii_unnormalized
    prefs.prenii_unnormalized_t1
    prefs.prenii_unnormalized_pd
    prefs.tranii_unnormalized
    prefs.cornii_unnormalized
    prefs.sagnii_unnormalized
    prefs.rawctnii_unnormalized
    prefs.ctnii_coregistered
    prefs.rest
    prefs.dti
};

extrafiles={
    prefs.prenii
    prefs.tranii
    prefs.cornii
    prefs.sagnii
    prefs.ctnii
    prefs.gprenii
    prefs.gtranii
    prefs.gcornii
    prefs.gsagnii
    prefs.gctnii
    prefs.tp_ctnii_coregistered
    prefs.tp_gctnii
    prefs.b0
    prefs.fa
    prefs.fa2anat
    ['l', prefs.fa2anat]
    ['gl', prefs.fa2anat]
    'grid.nii'
    'glgrid.nii'
};

leadfiles = [corefiles; extrafiles];
leadfiles = cellfun(@ea_niifileparts, leadfiles,'UniformOutput',0); % remove ".nii" extension
leadfiles = cellfun(@(x) x(3:end), leadfiles,'UniformOutput',0);    % remove "./" prefix from 'ea_niifileparts'
