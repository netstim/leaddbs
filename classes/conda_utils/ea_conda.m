classdef (Abstract) ea_conda

    properties (Constant)
        install_path = fullfile(ea_getearoot, 'ext_libs', 'miniconda');
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
            out = fullfile(ea_conda.install_path, bin_folder, ['conda' ext]);
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
                os_name = 'MacOSX';
            elseif ispc
                os_name = 'Windows';
            else
                os_name = 'Linux';
            end

            miniconda_installer_url =  ['https://repo.anaconda.com/miniconda/Miniconda3-latest-' os_name '-x86_64'];
            installer_file = fullfile(ea_conda.install_path, 'miniconda');

            if isunix
                miniconda_installer_url = [miniconda_installer_url '.sh'];
                installer_file = [installer_file '.sh'];
                install_call = ['bash ' installer_file ' -b -f -p ' ea_conda.install_path];
            else
                miniconda_installer_url = [miniconda_installer_url '.exe'];
                installer_file = [installer_file '.exe'];
                install_call = ['start /wait "" ' installer_file ' /InstallationType=JustMe /RegisterPython=0 /S /D=' ea_conda.install_path];
            end

            ea_conda.websave_verbose(installer_file, miniconda_installer_url);
            ea_conda.run_install_call(install_call)

            delete(installer_file);
            disp('Miniconda installed')

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
            disp('Installing miniconda...')
            [status,~] = system(install_call);
            if status
                error('Failed to install miniconda');
            end
        end

    end


end
