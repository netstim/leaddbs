function settings = ea_load_prob_parameter(options, settings, outputPaths, sideCode, sample_i)

% Load probabilistic sample
% IMPORTANT:this is a special case for multisource prob. PAM
% By Butenko, konstantinmgtu@gmail.com

arguments
    options             % Lead-DBS options for electrode reconstruction and stimulation
    settings            % parameters for OSS-DBS simulation
    outputPaths         % various paths to conform with lead-dbs BIDS structure 
    sideCode            % hemisphere suffix, 'rh' or 'lh' 
    sample_i            % index of the probabilistic sample
end

vatsettings = options.prefs.machine.vatsettings;
pparam = vatsettings.butenko_probabilistic_parameter;

load([outputPaths.outputDir, filesep, pparam,'_samples_',sideCode,'.mat'], 'prob_samples')
if strcmp(pparam,'Fiber Diameter')
    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
    settings.fiberDiameter(:) =  prob_samples(sample_i,1);
elseif strcmp(pparam,'Axon Length')
    settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
    settings.axonLength(:) = prob_samples(sample_i,1);
end
parameterFile = fullfile(outputPaths.outputDir, 'oss-dbs_parameters.mat');
save(parameterFile, 'settings', '-v7.3');