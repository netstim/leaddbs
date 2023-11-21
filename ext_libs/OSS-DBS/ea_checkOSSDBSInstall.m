function ea_checkOSSDBSInstall
% Check if OSS-DBS dependencies have been properly configured

if ~isempty(getenv('SINGULARITY_NAME')) % Singularity
    % Set pythonPath
    pythonPath = ea_findBinPath('python3');
    ea_setprefs('env.pythonPath', pythonPath, 'user');

    % Set installed flag
    prefs = ea_prefs;
    vatsettings = prefs.machine.vatsettings;
    vatsettings.oss_dbs.installed = 1;
    ea_setprefs('vatsettings', vatsettings);

    fprintf('OSS-DBS dependencies have been properly configured.\n');
    return
end

binPath = getenv('PATH'); % Current PATH

% Check docker image
ea_checkDocker('ningfei/oss-dbs:latest');

if isunix && ~ismac % Use local built image for Linux to solve permission issue
    [~, id] = system('docker images -q ningfei/oss-dbs:custom');
    if ~isempty(id)
        fprintf('docker image found: ningfei/oss-dbs:custom\n');
        fprintf('\nRebuilding local docker image...\n');
    else
        fprintf('\nBuilding local docker image...\n');
    end
    currentPath = pwd;
    cd([ea_getearoot, 'ext_libs/OSS-DBS']);
    system('docker build --build-arg UID=$(id -u) --build-arg GID=$(id -g) -t ningfei/oss-dbs:custom -f custom.Dockerfile .');
    cd(currentPath);
end

% Check OSS-DBS_Dependency installation
condaenv = ea_conda_env('OSS-DBS');
if ~condaenv.is_created
    ea_cprintf('CmdWinWarnings', 'Initializing OSS-DBS conda environment...\n')
    condaenv.create;
    ea_cprintf('CmdWinWarnings', 'OSS-DBS conda environment initialized.\n')
elseif ~condaenv.is_up_to_date
    ea_cprintf('CmdWinWarnings', 'Updating OSS-DBS conda environment...\n')
    condaenv.update;
    ea_cprintf('CmdWinWarnings', 'OSS-DBS conda environment initialized.\n')
end

fprintf('\n\nOSS-DBS dependencies installed under: %s\n', condaenv.path);

% Set PATH
if ispc
    pythonPath = condaenv.path;
    setenv('PATH', [pythonPath, ';', binPath]);
else
    pythonPath = fullfile(condaenv.path, 'bin');
    setenv('PATH', [pythonPath, ':', binPath]);
end

% Save prefs
ea_setprefs('env.pythonPath', pythonPath, 'user');

% Set installed flag
prefs = ea_prefs;
vatsettings = prefs.machine.vatsettings;
vatsettings.oss_dbs.installed = 1;
ea_setprefs('vatsettings', vatsettings);

fprintf('\nOSS-DBS dependencies have been properly configured.\n');
