function settings = ea_updatePAM_parameter(options,settings,outputPaths,sample_i)
% Update one PAM parameter that is subjected to uncertainty (set in the GUI)
% By Butenko, Rajamani and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    settings    % parameters for OSS-DBS simulation
    outputPaths % various paths to conform with lead-dbs BIDS structure
    sample_i  {mustBeNumeric} % sample index
end

% load the prob PAM parameters from the GUI
vatsettings = options.prefs.machine.vatsettings;
parameter_limits = vatsettings.butenko_parameter_limits; %parameter_limits = [2.0,4.0];  % interval
N_samples = vatsettings.butenko_N_samples;
parameter_step = (parameter_limits(2)-parameter_limits(1)) / (N_samples - 1);
sampling_distribution = vatsettings.butenko_sampling_distribution;

% get a sample
if strcmp(sampling_distribution,'Equidistant')
    param_to_change = parameter_limits(1) + (sample_i-1)*parameter_step;
elseif strcmp(sampling_distribution,'Uniform')
    % alternatively, it can be sampled uniformly
    param_to_change = parameter_limits(1) + rand * (parameter_limits(2)-parameter_limits(1));
elseif strcmp(sampling_distribution,'Normal')
    % or from a normal distribution
    mu = mean(parameter_limits);
    sigma = (mean(parameter_limits) - parameter_limits(1)) / 3;
    param_to_change = normrnd(mu,sigma);
end

% update the parameter
pparam = vatsettings.butenko_probabilistic_parameter;
if strcmp(pparam,'Fiber Diameter')
    settings.fiberDiameter = options.prefs.machine.vatsettings.butenko_fiberDiameter;
    settings.fiberDiameter(:) =  param_to_change;
elseif strcmp(pparam,'Axon Length')
    settings.axonLength = options.prefs.machine.vatsettings.butenko_axonLength;
    settings.axonLength(:) = param_to_change;
elseif strcmp(pparam,'Encapsulation Thickness')
    ea_warndlg("TBA!")
end

parameterFile = fullfile(outputPaths.outputDir, 'oss-dbs_parameters.mat');
save(parameterFile, 'settings', '-v7.3');
