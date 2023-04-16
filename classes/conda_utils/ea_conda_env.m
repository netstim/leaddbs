classdef ea_conda_env

    properties
        yml;
        name;
        path;
    end

    properties (Dependent)
        is_created;
    end

    properties (Access = private, Constant)
        conda_path = ea_conda.bin_file_path;
    end

    methods

        function obj = ea_conda_env(ymlname)
            obj.yml = fullfile(fileparts(mfilename('fullpath')), 'environments', [ymlname '.yml']);
            if ~isfile(obj.yml)
                ea_cprintf('CmdWinErrors', 'Environment yml file doesn''t exist!\n');
                obj.yml = [];
                return;
            end
            obj.name = regexp(fgetl(fopen(obj.yml)), '(?<=name:\s+).*', 'match', 'once');
            obj.path = fullfile(ea_conda.install_path, 'envs', obj.name);
        end

        function b = get.is_created(obj)
            b = isfolder(obj.path);
        end

        function force_create(obj)
            obj.remove;
            obj.create;
        end

        function remove(obj)
            system([ea_conda_env.conda_path ' env remove --name ' obj.name]);
        end

        function create(obj)
            disp(['Creating environment ' obj.name '...'])
            [status, cmdout] = system([ea_conda_env.conda_path ' env create -f ' obj.yml]);
            if status
                fprintf('%s\n', strtrim(cmdout));
                ea_cprintf('CmdWinErrors', 'Failed to create environment %s! Please check the log above.\n', obj.name)
            end
        end

        function run_script(obj, script_path)
            obj.system(['python ' script_path])
        end

        function system(obj, command)
            if ~obj.is_created
                error(['Create python environment ' obj.name ' from Lead-DBS menu to use this function']);
            end
            if isunix
                setup_command = ['export PATH=' fullfile(obj.path,'bin') ':$PATH;'];
            else
                setup_command = [fullfile(ea_conda.install_path, 'condabin', 'activate.bat') ' ' obj.name ' & '];
                command = obj.inject_exe_to_command(command);
            end
            system([setup_command command]);
        end
    end

    methods (Static, Access = private)

        function command = inject_exe_to_command(command)
            first_space = regexp(command,' ','once');
            if isempty(first_space)
                first_space = length(command)+1;
            end
            command = [command(1:first_space-1) '.exe' command(first_space:end)];
        end
    end
end

