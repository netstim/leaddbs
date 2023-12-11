function ea_checkOSSDBSInstallv2

% OSS-DBS v2 deployment

% Check conda installation
if ~ea_conda.is_installed
    ea_conda.install;
end

% optional: provide the local path to OSS-DBSv2 in OSS-DBS-v2.yml if rep is not avaialble

% install OSS-DBS v2 in the virtual environment
env = ea_conda_env('OSS-DBS-v2.yml');
env.create;

if ~isunix
    ea_warndlg("To run PAM on Windows, please install NEURON from https://www.neuron.yale.edu/neuron/download")
end

% set installed flag
prefs = ea_prefs;
vatsettings = prefs.machine.vatsettings;
vatsettings.oss_dbs.installed = 1;
ea_setprefs('vatsettings', vatsettings);