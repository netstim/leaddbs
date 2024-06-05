function settings = ea_prepare_ossdbs(options)
% Import settings from Lead-DBS.
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options % Lead-DBS options for electrode reconstruction and stimulation
end

% Check OSS-DBS installation, set env
env = ea_conda_env('OSS-DBSv2');
ea_checkOSSDBSInstallv2(env);

% Set python path
binPath = getenv('PATH');
if isunix
    pythonPath = [env.path, filesep, 'bin'];
    setenv('PATH', [pythonPath, ':', binPath]);
else
    pythonPath = [env.path,';',env.path,filesep,'Scripts'];
    setenv('PATH', [pythonPath, ';', binPath]);
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

%% Lead-DBS hardwired parameters

% Set patient folder
settings.Patient_folder = options.subj.subjDir;
settings.Electrode_type = options.elmodel;
