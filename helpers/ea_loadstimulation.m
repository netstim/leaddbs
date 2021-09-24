function S=ea_loadstimulation(stimname,options)

stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimname);
filePrefix = ['sub-', options.subj.subjId, '_desc-'];

load([stimDir, filesep, filePrefix, 'stimparameters.mat']);
