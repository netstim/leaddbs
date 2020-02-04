function S=ea_loadstimulation(stimname,options)

load([options.root,options.patientname,filesep,'stimulations',filesep,ea_nt(options),stimname,filesep,'stimparameters.mat']);
