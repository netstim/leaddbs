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

% Check docker installation
dockerPath = ea_findBinPath('docker');
if isempty(dockerPath)
    if ismac
        % /usr/local/bin might not be in the system path env on macOS
        ea_error(sprintf(['docker not found!\nIf it''s already installed, ', ...
            'please run the line below in your terminal (not MATLAB terminal!) and reboot:\n', ...
            'sudo launchctl config user path /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin']), ...
            'Error', dbstack, 0);
    else
        ea_error('docker not found!', showdlg = 0, simpleStack = 1);
    end
else
    if ismac || ispc % Use upstream image since there's no permission issue
        [~, id] = system('docker images -q ningfei/oss-dbs:latest');
        if ~isempty(id)
            fprintf('docker image found: ningfei/oss-dbs:latest\n');
        end
        fprintf('\nPulling docker image...\n'); % Always pull to update local image
        system('docker pull ningfei/oss-dbs:latest');
    else % Use local built image
        [~, id] = system('docker images -q ningfei/oss-dbs:custom');
        if ~isempty(id)
            fprintf('docker image found: ningfei/oss-dbs:custom\n');
            fprintf('\nRebuilding docker image...\n');
        else
            fprintf('\nBuilding docker image...\n');
        end
        currentPath = pwd;
        cd([ea_getearoot, 'ext_libs/OSS-DBS']);
        system('docker pull ningfei/oss-dbs:latest'); % Pull to update base image
        system('docker build --build-arg UID=$(id -u) --build-arg GID=$(id -g) -t ningfei/oss-dbs:custom -f custom.Dockerfile .');
        cd(currentPath);
    end
end

% Check conda installation
if ~ea_conda.is_installed
    ea_conda.install;
end

% Check OSS-DBS_Dependency installation
condaenv = ea_conda_env('OSS-DBS');
if ~condaenv.is_created
    condaenv.create;
end

fprintf('\nOSS-DBS dependencies installed under: %s\n', condaenv.path);

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
