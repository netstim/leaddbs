function [settings,S] = ea_prepare_ossdbs(options,S)
% Import settings from Lead-DBS.
% For StimSets also add a dummy S
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options % Lead-DBS options for electrode reconstruction and stimulation
    S       % lead-dbs stimulation settings
end

% Check OSS-DBS installation, set env
if (~isfield(options.prefs,'ext_oss_env') || strcmp(options.prefs.ext_oss_env,'None')) && (~isfield(options.prefs,'wsl_env') || ~options.prefs.wsl_env)
    env = ea_conda_env('OSS-DBSv2');
    ea_checkOSSDBSInstallv2(env);
end

% Double check if lead is supported by OSS-DBS.
if ~ismember(options.elmodel, ea_ossdbs_elmodel)
    ea_error([options.elmodel, 'is not supported by OSS-DBS yet!'], simpleStack = 1);
end

%% settings GUI parameters
settings.butenko_segmAlg = options.prefs.machine.vatsettings.butenko_segmAlg;
settings.butenko_intersectStatus = options.prefs.machine.vatsettings.butenko_intersectStatus;
settings.removeElectrode = options.prefs.machine.vatsettings.butenko_removeElectrode;
settings.neuronModel = options.prefs.machine.vatsettings.butenko_axonModel;
settings.signalType = options.prefs.machine.vatsettings.butenko_signalType;
settings.biphasic = options.prefs.machine.vatsettings.butenko_biphasic;
settings.butenko_tensorData = options.prefs.machine.vatsettings.butenko_tensorData;
settings.AdaptiveRef = options.prefs.machine.vatsettings.butenko_AdaptiveRef;
settings.encapsulationType = options.prefs.machine.vatsettings.butenko_encapsulation;
settings.adaptive_threshold = options.prefs.machine.vatsettings.butenko_adaptive_ethresh;
% check what we simulate
settings.calcAxonActivation = options.prefs.machine.vatsettings.butenko_calcPAM;

if ~isfield(options.prefs.machine.vatsettings,'butenko_cond_model')
    settings.cond_model = 'ColeCole4';
elseif strcmp(options.prefs.machine.vatsettings.butenko_cond_model,'Homogeneous')
    settings.cond_model = 'Constant';
else
    settings.cond_model = options.prefs.machine.vatsettings.butenko_cond_model;
end

try
    % compute PAM over an uncertain parameter (set in the GUI)
    settings.prob_PAM = options.prefs.machine.vatsettings.butenko_prob_PAM;
    if settings.prob_PAM
        settings.N_samples = options.prefs.machine.vatsettings.butenko_N_samples;
    else
        settings.N_samples = 1;
    end
catch
    settings.prob_PAM = 0;
    settings.N_samples = 1;
end

settings.exportVAT = options.prefs.machine.vatsettings.butenko_calcVAT;
% Set native/MNI flag
settings.Estimate_In_Template = ~options.native;
% Set stimSetMode flag
settings.stimSetMode = options.stimSetMode;

% advance options
if isfield(options,'optimizer')
    % check if OSS-DBS is called via ea_OSS_optimizer
    settings.optimizer = options.optimizer;
    if settings.optimizer
        S_true_label = S.label;
        S = ea_set_optimizer(options,[options.subj.stimDir, filesep, ea_nt(options.native), S.label]);
        S.label = S_true_label;
        if settings.exportVAT
            % special case of VTA-based optimization
            options.PathwayTune_master_dict = [];
        elseif settings.calcAxonActivation
            settings.PathwayTune_master_dict = options.PathwayTune_master_dict;
        end

    end
else
    settings.optimizer = 0;
end

if settings.stimSetMode == 1 && ~settings.optimizer
    % default stimulation (Stim_protocols will be actually used)
    S.Rs1.amp = 1.0;
    S.Ls1.amp = 1.0;
    S.amplitude{1,1}(1) = 1.0;
    S.amplitude{1,2}(1) = 1.0;
end

try
    % check if OSS-DBS is called for ANN training (StimSet case)
    settings.trainANN = options.trainANN;
catch
    settings.trainANN = 0;
end

if settings.optimizer || settings.trainANN
    settings.stimSetMode = 1; % OSS-DBS will solve a "unit" problem
end

if (settings.optimizer || settings.trainANN || settings.stimSetMode) && strcmp(settings.butenko_intersectStatus,'activated_at_active_contacts')
    ea_error("Option 'Activated near active contacts' is not supported for StimSet mode")
end

if settings.stimSetMode
    settings.current_control = [1;1];
end

%% Lead-DBS hardwired parameters

% Set patient folder
settings.Patient_folder = options.subj.subjDir;
settings.Electrode_type = options.elmodel;
