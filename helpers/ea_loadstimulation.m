function S=ea_loadstimulation(stimname,options)

load([options.root,options.patientname,filesep,'stimulations',filesep,stimname,filesep,'stimparameters.mat']);
