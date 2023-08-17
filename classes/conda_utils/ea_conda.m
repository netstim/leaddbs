classdef (Abstract) ea_conda

    properties (Constant)
        install_path = fullfile(ea_getearoot, 'ext_libs', 'mambaforge');
    end

    methods (Static)

        function out = bin_file_path
            if isunix
                bin_folder = 'bin';
                ext = '';
            else
                bin_folder = 'condabin';
                ext = '.bat';
            end
            out = fullfile(ea_conda.install_path, bin_folder, ['mamba' ext]);
        end

        function b = is_installed
            b = isfile(ea_conda.bin_file_path);
        end

        function remove
            rmdir(ea_conda.install_path, 's')
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

            mambaforge_installer_url =  ['https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-' osarch];
            installer_file = fullfile(ea_getearoot, 'mambaforge');

            if isunix
                mambaforge_installer_url = [mambaforge_installer_url '.sh'];
                installer_file = [installer_file '.sh'];
                install_call = ['bash ' installer_file ' -b -f -p ' ea_conda.install_path];
            else
                mambaforge_installer_url = [mambaforge_installer_url '.exe'];
                installer_file = [installer_file '.exe'];
                install_call = ['start /wait "" ' installer_file ' /InstallationType=JustMe /RegisterPython=0 /S /D=' ea_conda.install_path];
            end

            ea_conda.websave_verbose(installer_file, mambaforge_installer_url);
            ea_conda.run_install_call(install_call)

            delete(installer_file);
            disp('mambaforge installed')
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

        function update_base
            conda = ea_conda.bin_file_path;
            system([conda, ' update conda mamba -y']);
            system([conda, ' update --all -y']);
            fprintf('\n');
        end
    end

    methods (Access = private, Static)

        function websave_verbose(filename, url)
            disp(['Downloading ' url ' to ' filename]);
            try
                websave(filename, url);
            catch
                error(['Failed to download ' url]);
            end
        end

        function run_install_call(install_call)
            disp('Installing mambaforge...')
            [status,~] = system(install_call);
            if status
                error('Failed to install mambaforge');
            end
        end
    end
end
