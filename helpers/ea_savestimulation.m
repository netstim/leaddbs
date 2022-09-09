function ea_savestimulation(S,options)

stimDir = fullfile(options.subj.stimDir, ea_nt(options), S.label);
ea_mkdir(stimDir);

filePrefix = ['sub-', options.subj.subjId, '_desc-'];

save([stimDir, filesep, filePrefix, 'stimparameters.mat'], 'S');
