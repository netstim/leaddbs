function settings = ea_prepare_ossdbs(options)

% Check OSS-DBS installation, set env
env = ea_conda_env('OSS-DBSv2');
ea_checkOSSDBSInstallv2(env);

binPath = getenv('PATH'); % Current PATH
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
settings.outOfCore = 0;
% check what we simulate
settings.calcAxonActivation = options.prefs.machine.vatsettings.butenko_calcPAM;
settings.exportVAT = options.prefs.machine.vatsettings.butenko_calcVAT;

%% Lead-DBS hardwired parameters

% Set patient folder
settings.Patient_folder = options.subj.subjDir;

% Set native/MNI flag
settings.Estimate_In_Template = ~options.native;

settings.Electrode_type = options.elmodel;