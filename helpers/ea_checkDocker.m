function ea_checkDocker(image)
% Function to check docker installation [and docker image].
% If docker is installed, the specified image will be pulled.

arguments
    image {mustBeTextScalar} = ''
end

% First check docker installation
dockerPath = ea_findBinPath('docker');
if isempty(dockerPath)
    if ismac
        % /usr/local/bin might not be in the system path env on macOS
        ea_error(sprintf(['docker not found!\nIf it''s already installed, ', ...
            'please run the line below in your terminal (NOT MATLAB Command Window!) and reboot:\n', ...
            'sudo launchctl config user path /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin']), ...
            showdlg = 0, simpleStack = 1);
    else
        ea_error('docker not found!', showdlg = 0, simpleStack = 1);
    end
else
    fprintf('docker found: %sdocker\n', [dockerPath, filesep]);
end

% Check docker image
if ~isempty(image)
    [~, id] = system(['docker images -q ' image]);
    if ~isempty(id)
        fprintf('docker image found: %s\n', image);
    end
    fprintf('\nPulling docker image...\n'); % Always pull to update local image
    system(['docker pull ' image]);
end
