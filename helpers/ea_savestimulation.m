function ea_savestimulation(S,options)

if ~isfield(options,'root') % called from lead group
    return
end
if ~exist([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),S.label],'dir')
    mkdir([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),S.label]);
end

save([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),S.label,filesep,'stimparameters.mat'],'S');
