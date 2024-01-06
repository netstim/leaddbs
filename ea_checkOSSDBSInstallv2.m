function ea_checkOSSDBSInstallv2

% OSS-DBS v2 deployment

% Check conda installation
if ~ea_conda.is_installed
    ea_cprintf('*Comments', 'Initializing conda...\n');
    ea_conda.install;
    ea_cprintf('*Comments', 'Done...\n');
end

% optional: provide the local path to OSS-DBSv2 in OSS-DBSv2.yml if rep is not avaialble

% install OSS-DBS v2 in the virtual environment
env = ea_conda_env('OSS-DBSv2');
if ~env.is_up_to_date
    ea_cprintf('*Comments', 'Updating OSS-DBS v2 conda environment...\n');
    env.force_create;
    ea_cprintf('*Comments', 'Done.\n');
end

if ~isunix
    msgbox(sprintf('An external installer for NEURON will be opened.\nPlease install it using the default parameters.'), '', 'help', 'modal');
    installer = fullfile(ea_prefsdir, 'temp', 'nrn-8.2.3.exe');
    ea_mkdir(fileparts(installer));
    try
        websave(installer, 'https://github.com/neuronsimulator/nrn/releases/download/8.2.3/nrn-8.2.3.w64-mingw-py-37-38-39-310-311-setup.exe');
    catch ME
        ea_error(['Failed to download NEURON installer for Windows:\n', ME.message], simpleStack=true);
    end
    system(installer);
    ea_delete(installer);
else
    env.system('pip3 install neuron==8.2.3')
end

% set installed flag
prefs = ea_prefs;
vatsettings = prefs.machine.vatsettings;
vatsettings.oss_dbs.installed = 1;
ea_setprefs('vatsettings', vatsettings);