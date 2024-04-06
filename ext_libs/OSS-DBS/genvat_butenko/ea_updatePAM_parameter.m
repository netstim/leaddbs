function settings = ea_updatePAM_parameter(options,settings,N_samples,outputPaths,sample_i)
% Update one PAM parameter that is subjected to uncertainty, for example,
% fiber diameter
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    settings    % parameters for OSS-DBS simulation
    N_samples {mustBeNumeric} % total number of samples
    outputPaths % various paths to conform with lead-dbs BIDS structure
    sample_i  {mustBeNumeric} % sample index
end

% hardcoded for now
parameter_limits = [2.0,4.0];
parameter_step = (parameter_limits(2)-parameter_limits(1)) / (N_samples - 1);

% alternatively, can be sampled randomly
settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
settings.fiberDiameter(:) = parameter_limits(1);
settings.fiberDiameter(:) = settings.fiberDiameter(:) + (sample_i-1)*parameter_step;

parameterFile = fullfile(outputPaths.outputDir, 'oss-dbs_parameters.mat');
save(parameterFile, 'settings', '-v7.3');