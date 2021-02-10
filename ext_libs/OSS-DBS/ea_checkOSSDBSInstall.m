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
            'please run the line below in your terminal and reboot:\n', ...
            'sudo launchctl config user path /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin']), ...
            'Error', dbstack, 0);
    else
        ea_error('docker not found!', 'Error', dbstack, 0);
    end
else
    if ismac || ispc % Use upstream image since there's no permission issue
        [~, id] = system('docker images -q ningfei/oss-dbs');
        if ~isempty(id)
            fprintf('docker image found: ningfei/oss-dbs\n');
        end
        fprintf('\nPulling docker image...\n'); % Always pull to update local image
        system('docker pull ningfei/oss-dbs:latest');
    else % Use local built image
        [~, id] = system('docker images -q custom_oss-dbs');
        if ~isempty(id)
            fprintf('docker image found: custom_oss-dbs\n');
            fprintf('\nRebuilding docker image...\n');
        else
            fprintf('\nBuilding docker image...\n');
        end
        currentPath = pwd;
        cd([ea_getearoot, 'ext_libs/OSS-DBS']);
        system('docker pull ningfei/oss-dbs:latest'); % Pull to update base image
        system('docker build --build-arg UID=$(id -u) --build-arg GID=$(id -g) -t custom_oss-dbs .');
        cd(currentPath);
    end
end

% Check if there's user defined python path
prefs = ea_prefs;
if isfield(prefs.env, 'pythonPath')
    pythonPath = prefs.env.pythonPath;
    if isunix
        setenv('PATH', [pythonPath, ':', binPath]);
    else
        setenv('PATH', [pythonPath, ';', binPath]);
    end
end

% Check python3
if ispc
    pythonBinName = 'python';
else
    pythonBinName = 'python3';
end
pythonPath = ea_findBinPath(pythonBinName);

fprintf('\n');
if isempty(pythonPath)
    ea_error('python3 not found!', 'Error', dbstack, 0);
else
    % Confirm python3 path
    answer = questdlg(['python3 found under ', pythonPath],...
        'python3 found...',...
        'Confirm',...
        'Customize',...
        'Confirm');

    % Use custom python3
    if strcmp(answer, 'Customize')
        if isunix
            pythonPath = ea_uigetdir('/usr/local/bin', 'Select python3 path...');
        else
            pythonPath = ea_uigetdir('C:\Program Files', 'Select python3 path...');
        end
        pythonPath = pythonPath{1};
    end
end

% python3 executable
pythonBinPath = [pythonPath, filesep, pythonBinName];
fprintf('python3 detected: %s\n', pythonBinPath);

% Check h5py
[status, h5pyPath] = system([pythonBinPath, ' -c "import h5py;print(h5py.__file__)"']);
if status
    ea_error(['h5py not found! Please run ''', ...
             pythonBinPath,' -m pip install h5py',...
             ''' in your terminal.'], 'Error', dbstack, 0);
else
    fprintf('h5py detected: %s\n', fileparts(h5pyPath));
end

% Check PyQt5
[status, pyqt5Path] = system([pythonBinPath, ' -c "import PyQt5;print(PyQt5.__file__)"']);
if status
    ea_error(['PyQt5 not found! Please run ''', ...
             pythonBinPath,' -m pip install PyQt5',...
             ''' in your terminal.'], 'Error', dbstack, 0);
else
    fprintf('PyQt5 detected: %s\n', fileparts(pyqt5Path));
end

% All passed, set env, save prefs
if isunix
    setenv('PATH', [pythonPath, ':', binPath]);
else
    setenv('PATH', [pythonPath, ';', binPath]);
end
ea_setprefs('env.pythonPath', pythonPath, 'user');

% Set installed flag
vatsettings = prefs.machine.vatsettings;
vatsettings.oss_dbs.installed = 1;
ea_setprefs('vatsettings', vatsettings);

fprintf('\nOSS-DBS dependencies have been properly configured.\n');
