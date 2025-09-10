classdef ea_ext_env

    properties
        name;
        path;
        python;
    end

    properties (Dependent)
        is_created;
        installed_version;
        latest_version;
    end

    methods

        function obj = ea_ext_env(env_path)

            [path_to_env,obj.name,~] = fileparts(env_path);
            obj.path = env_path;
            if isunix
                obj.python = fullfile(env_path, 'bin', 'python');
            else
                obj.python = fullfile(env_path, 'python.exe');
            end
        end

        function b = get.is_created(obj)
            b = isfolder(obj.path);
        end

        function run_script(obj, script_path)
            obj.system(['python ' script_path])
        end

        function varargout = system(obj, command)
            if ~obj.is_created
                error(['Create external python environment ' obj.name]);
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
