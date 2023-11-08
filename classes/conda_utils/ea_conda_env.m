classdef ea_conda_env

    properties
        name;
        yml;
        path;
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
            ymlfile = ea_regexpdir(fullfile(fileparts(mfilename('fullpath')), 'environments'), ymlname, 0);
            if ~isempty(ymlfile)
                if length(ymlfile) > 1
                    ea_cprintf('CmdWinWarnings', ...
                        'Duplicated environment files found. Consider cleaning up the folder:\n%s\n', ...
                        fullfile(fileparts(mfilename('fullpath')), 'environments'));
                end
                obj.yml = ymlfile{1};
            else
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

        function ver = get.installed_version(obj)
            verFile = ea_regexpdir(obj.path, '^v\d+', 0, 'f');
            if isempty(verFile)
                ver = '';
            else
                [~, verFile] = fileparts(verFile{1});
                ver = verFile(2:end);
            end
        end

        function ver = get.latest_version(obj)
            yaml = readyaml(obj.yml);
            if ~isfield(yaml, 'version')
                ver = '';
                ea_cprintf('CmdWinWarnings', 'Missing version tag in env yaml definition.\n');
            else
                ver = num2str(yaml.version);
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
            system([obj.mamba_path ' env remove --name ' obj.name]);
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
            [status, cmdout] = system([obj.mamba_path ' env create -f ' obj.yml]);
            if status
                fprintf('%s\n', strtrim(cmdout));
                ea_cprintf('CmdWinErrors', 'Failed to create environment %s! Please check the log above.\n', obj.name)
            end
        end

        function run_script(obj, script_path)
            obj.system(['python ' script_path])
        end

        function status = system(obj, command)
            if ~obj.is_created
                error(['Create python environment ' obj.name ' from Lead-DBS menu to use this function']);
            end
            if isunix
                setup_command = ['export PATH=' fullfile(obj.path,'bin') ':$PATH;'];
            else
                setup_command = [fullfile(ea_conda.install_path, 'condabin', 'activate.bat') ' ' obj.name ' & '];
                command = obj.inject_exe_to_command(command);
            end
            status = system([setup_command command]);
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

