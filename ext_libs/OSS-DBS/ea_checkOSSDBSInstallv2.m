function ea_checkOSSDBSInstallv2(env)

% Check conda installation
if ~ea_conda.is_installed
    ea_cprintf('*Comments', 'Initializing conda...\n');
    ea_conda.install;
    ea_cprintf('*Comments', 'Done...\n');
end

if ~exist('env', 'var')
    env = ea_conda_env('OSS-DBSv2');
end

% Install/Update OSS-DBS v2 conda environment
if ~env.is_up_to_date
    ea_cprintf('*Comments', 'Updating OSS-DBS v2 conda environment...\n');
    env.force_create;
    ea_cprintf('*Comments', 'Done.\n');
end

% Check NEURON install
version = '8.2.4';
if ispc
    [status, cmdout] = system('neuron --version');
    if status || ~contains(cmdout, version)
        ea_cprintf('*Comments', 'Installing NEURON %s for Windows...\n', version);
        installer = fullfile(ea_prefsdir, 'temp', ['nrn-', version, '.exe']);
        ea_mkdir(fileparts(installer));
        try
            websave(installer, ['https://github.com/neuronsimulator/nrn/releases/download/', version, '/nrn-', version, '.w64-mingw-py-37-38-39-310-311-setup.exe']);
        catch ME
            ea_error(['Failed to download NEURON installer for Windows:\n', ME.message], simpleStack=true);
        end
        installFolder = fullfile(ea_prefsdir, 'neuron');
        ea_delete(installFolder);
        try
            system(['start /b /wait "Install NEURON" ', path_helper(installer), ' /S /D=', path_helper(installFolder)]);
        catch ME
            ea_error(['Failed to install NEURON for Windows:\n', ME.message], simpleStack=true);
        end
        ea_delete(installer);
        setenv('PATH', [getenv('PATH'), ';', fullfile(installFolder, 'bin')]);
    end
else
    [status, cmdout] = env.system('python -c ''import neuron;print(neuron.__version__)''');
    if status || ~contains(cmdout, version)
        ea_cprintf('*Comments', 'Install NEURON %s...\n', version);
        env.system(['pip3 install -U neuron==', version]);
    end
end


% Handle space in path on Windows
function path = path_helper(path)
parts = strsplit(path, filesep);
for i=1:length(parts)
    if contains(parts{i}, ' ')
        parts{i} = ['"' parts{i} '"'];
    end
end
path = strjoin(parts, filesep);

