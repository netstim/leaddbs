classdef (Abstract) ea_conda

    methods (Static)
        function path = install_path(custompath)
            if exist('custompath', 'var') && ~isempty(custompath)
                % custom installation path
                path = custompath;
                ea_setprefs('conda.install_path', path, 'user');
                return;
            end
            prefs = ea_prefs;
            if ~isempty(prefs.conda.install_path)
                % Installation path already set in prefs
                path = prefs.conda.install_path;
            else
                % Set installation path
                path = fullfile(ea_prefsdir, 'miniforge');
                answer = questdlg(sprintf('Confirm Conda installation folder:\n%s', path), '', 'Yes', 'Custom', 'Yes');
                if strcmp(answer, 'Custom')
                    custompath = uigetdir(path, 'Specify Conda installation folder');
                    if ischar(custompath)
                        if ~endsWith(custompath, {'conda', 'condaforge', 'miniforge', 'mambaforge'}, 'IgnoreCase', true)
                            path = fullfile(custompath, 'miniforge');
                        else
                            path = custompath;
                        end
                    end
                end
                ea_setprefs('conda.install_path', path, 'user');
            end
        end

        function path = mamba_path
            if isunix
                bin_folder = 'bin';
                ext = '';
            else
                bin_folder = 'condabin';
                ext = '.bat';
            end
            path = fullfile(ea_conda.install_path, bin_folder, ['mamba' ext]);
        end

        function b = is_installed
            b = isfile(ea_conda.mamba_path);
        end

        function install
            mkdir(ea_conda.install_path);

            if ismac
                [~, arch] = system('uname -a');
                if contains(arch, 'arm64', 'IgnoreCase', true)
                    osarch = 'MacOSX-arm64';
                else
                    osarch = 'MacOSX-x86_64';
                end
            elseif ispc
                osarch = 'Windows-x86_64';
            else
                osarch = 'Linux-x86_64';
            end

            miniforge_installer_url =  ['https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-' osarch];
            installer_file = fullfile(ea_prefsdir, 'miniforge');

            if isunix
                miniforge_installer_url = [miniforge_installer_url '.sh'];
                installer_file = [installer_file '.sh'];
                install_call = ['bash ' installer_file ' -b -f -p ' ea_conda.install_path];
            else
                miniforge_installer_url = [miniforge_installer_url '.exe'];
                installer_file = [installer_file '.exe'];
                install_call = ['start /wait "" ' installer_file ' /InstallationType=JustMe /RegisterPython=0 /S /D=' ea_conda.install_path];
            end

            ea_conda.websave_verbose(installer_file, miniforge_installer_url);
            ea_conda.run_install_call(install_call)

            delete(installer_file);
            ea_cprintf('*Comments', 'miniforge installed...\n');

            % Set some conda configs
            [~, cmdout] = ea_conda.run('conda config --get channels');
            if isempty(cmdout)
                ea_conda.run('conda config --prepend channels conda-forge');
                ea_conda.run('conda config --append channels defaults');
            else
                if ~contains(cmdout, "channels 'conda-forge'")
                    ea_conda.run('conda config --prepend channels conda-forge');
                end
                if ~contains(cmdout, "channels 'defaults'")
                    ea_conda.run('conda config --append channels conda-forge');
                end
            end

            [~, cmdout] = ea_conda.run('conda config --get ssl_verify');
            if isempty(cmdout)
                ea_conda.run('conda config --set ssl_verify false');
            end

            [~, cmdout] = ea_conda.run('conda config --get auto_activate_base');
            if isempty(cmdout)
                ea_conda.run('conda config --set auto_activate_base false');
            end

            ea_cprintf('*Comments', 'Please run ''ea_conda_setproxy'' if you are behind a proxy.\n');
        end

        % Run command (for example mamba or pip) in conda base environment
        function varargout = run(command)
            if ispc
                bin_folder = fullfile(ea_conda.install_path, 'Scripts');
                setenv('PATH', [bin_folder ';' getenv('PATH')]);
            else
                bin_folder = fullfile(ea_conda.install_path, 'bin');
                setenv('PATH', [bin_folder ':' getenv('PATH')]);
            end

            if nargout <= 1
                varargout{1} = system(command);
            else
                [varargout{1}, varargout{2}]  = system(command);
                varargout{2} = strip(varargout{2});
            end
        end

        function update_base
            conda = ea_conda.mamba_path;
            system([conda, ' update conda mamba -y']);
            system([conda, ' update --all -y']);
            fprintf('\n');
        end

        function clean
            conda = ea_conda.mamba_path;
            system([conda, ' clean -tpy']);
            fprintf('\n');
        end

        function remove
            ea_delete(ea_conda.install_path)
        end

        function listenv
            ymlFolder = fullfile(ea_getearoot, 'classes', 'conda_utils', 'environments');
            ymlFile = ea_regexpdir(ymlFolder, '\.yml', 0, 'f');
            envName = regexp(ymlFile, ['(?<=environments\', filesep, ').*(?=\.yml$)'], 'match', 'once');
            for i=1:length(envName)
                env = ea_conda_env(envName{i});
                if env.is_created
                    envName{i} = [envName{i}, ' (Installed)'];
                end
            end
            fprintf('%s\n', strjoin(envName, '\n'));
        end
    end

    methods (Access = private, Static)
        function websave_verbose(filename, url)
            ea_cprintf('*Comments', 'Downloading %s to %s\n', url, filename);
            try
                websave(filename, url);
            catch
                error(['Failed to download ' url]);
            end
        end

        function run_install_call(install_call)
            ea_cprintf('*Comments', 'Installing miniforge...\n');
            [status,~] = system(install_call);
            if status
                error('Failed to install miniforge');
            end
        end
    end
end
