classdef ea_conda_env

    properties
        name;
        yml;
        path;
        python;
    end

    properties (Dependent)
        is_created;
        installed_version;
        latest_version;
    end

    properties (Access = private, Constant)
        mamba_path = ea_conda.mamba_path;
    end

    methods

        function obj = ea_conda_env(ymlname)
            ymlname = erase(ymlname, '.yml');
            ymlfile = ea_regexpdir(fullfile(fileparts(mfilename('fullpath')), 'environments'), ['^', ymlname, '\.yml$'], 0);
            if ~isempty(ymlfile)
                obj.yml = ymlfile{1};
            else
                ea_cprintf('CmdWinErrors', 'Environment yml file doesn''t exist!\n');
                obj.yml = [];
                return;
            end
            fileID = fopen(obj.yml);
            obj.name = regexp(fgetl(fileID), '(?<=name:\s+).*', 'match', 'once');
            fclose(fileID);
            obj.path = fullfile(ea_conda.install_path, 'envs', obj.name);
            if isunix
                obj.python = fullfile(obj.path, 'bin', 'python');
            else
                obj.python = fullfile(obj.path, 'python.exe');
            end
        end

        function b = get.is_created(obj)
            b = isfolder(obj.path);
        end

        function ver = get.installed_version(obj)
            try
                state = loadjson(fullfile(obj.path, 'conda-meta', 'state'));
                ver = num2str(state.env_vars.env_version);
            catch
                ver = '';
            end
        end

        function ver = get.latest_version(obj)
            yaml = readyaml(obj.yml);
            try
                ver = num2str(yaml.variables.env_version);
            catch
                ver = '';
                ea_cprintf('CmdWinWarnings', 'Missing version tag in env yaml definition.\n');
            end
        end

        function up_to_date = is_up_to_date(obj)
            up_to_date = strcmp(obj.installed_version, obj.latest_version);
        end

        function force_create(obj)
            obj.remove;
            obj.create;
        end

        function update(obj)
            obj.force_create;
        end

        function remove(obj)
            if isfile(obj.mamba_path)
                system([obj.mamba_path ' env remove --yes --name ' obj.name]);
            elseif isfolder(obj.path)
                ea_cprintf('CmdWinWarnings', 'Missing mamba binary! Deleting %s env folder directly.\n', obj.name);
            end
            ea_delete(obj.path);
        end

        function create(obj)
            if isfolder(obj.path)
                ea_cprintf('CmdWinErrors', 'Conda env installation folder already exists!\nConsider using ''update''.\n');
                return;
            end

            if ~isfile(obj.mamba_path)
                ea_cprintf('CmdWinWarnings', 'Conda installation not found! Installing now...\n');
                ea_conda.install;
            end

            disp(['Creating environment ' obj.name '...'])
            obj.yml
            
            [status, cmdout] = system([obj.mamba_path ' env create -f ' obj.yml]);
            if status
                fprintf('%s\n', strtrim(cmdout));
                ea_cprintf('CmdWinErrors', 'Failed to create environment %s! Please check the log above.\n', obj.name)
            else
                system([obj.mamba_path ' clean -tpyq']);
            end
        end

        function run_script(obj, script_path)
            obj.system(['python ' script_path])
        end

        function varargout = system(obj, command)
            if ~obj.is_created
                error(['Create python environment ' obj.name ' from Lead-DBS menu to use this function']);
            end

            pathEnv = getenv('PATH');

            if isunix
                setenv('PATH', [fullfile(obj.path, 'bin') ':' getenv('PATH')]);
            else
                setenv('PATH', [obj.path ';' fullfile(obj.path, 'Scripts') ';' getenv('PATH')]);
            end

            if nargout <= 1
                varargout{1} = system(command);
            else
                [varargout{1}, varargout{2}] = system(command);
                varargout{2} = strip(varargout{2});
            end

            % Restore env after running conda command
            setenv('PATH', pathEnv)
        end
    end
end
