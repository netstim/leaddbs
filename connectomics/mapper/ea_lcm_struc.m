function ea_lcm_struc(options)

if strcmp(options.lcm.struc.connectome,'No structural connectome found.')
    return
end
disp('Running structural connectivity...');

useNativeSeed = options.prefs.lcm.struc.patienttracts.nativeseed;
if strcmp(options.lcm.struc.connectome, 'Patient''s fiber tracts')
    connBaseFolder = 'Patient''s fiber tracts';
    if useNativeSeed
        options.lcm.struc.connectome = strrep(options.prefs.FTR_unnormalized, '.mat', '_anat.mat');
    else
        options.lcm.struc.connectome = options.prefs.FTR_normalized;
    end
else
    connBaseFolder = ea_getconnectomebase();
end

cs_dmri_conseed(connBaseFolder, options);
disp('Done.');
