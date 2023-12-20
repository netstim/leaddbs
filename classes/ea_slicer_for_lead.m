classdef ea_slicer_for_lead

    properties
        install_path
        executable_path
        installed_version
        upstream_version
        download_url
    end

    properties (Access = private)
        release_tag = 'v231122-beta'
        release_version = struct('linux', '0.1.0-2023-11-22', 'mac', '0.1.0-2023-11-22', 'win', '0.1.0-2023-11-22')
    end

    methods
        function obj = ea_slicer_for_lead()
            obj.install_path = fullfile(ea_prefsdir, 'SlicerForLeadDBS');

            if ismac
                obj.executable_path = fullfile(obj.install_path, 'SlicerForLeadDBS.app', 'Contents', 'MacOS', 'SlicerForLeadDBS');
            elseif ispc
                obj.executable_path = fullfile(obj.install_path, 'SlicerForLeadDBS.exe');
            else
                obj.executable_path = fullfile(obj.install_path, 'SlicerForLeadDBS');
            end

            ver_file = ea_regexpdir(obj.install_path, '^v.*', 0, 'f');
            if ~isempty(ver_file)
                obj.installed_version = erase(ver_file{1}, [obj.install_path, filesep, 'v']);
            else
                obj.installed_version = '';
            end

            if ismac
                os = 'macosx';
                arch = 'amd64';
                obj.upstream_version = obj.release_version.mac;
            elseif isunix
                os = 'linux';
                arch = 'amd64';
                obj.upstream_version = obj.release_version.linux;
            elseif ispc
                os = 'win';
                arch = 'amd64';
                obj.upstream_version = obj.release_version.win;
            end
            download_file = ['SlicerForLeadDBS-' obj.upstream_version '-' os '-' arch '.zip'];
            obj.download_url = ['https://github.com/netstim/SlicerForLeadDBS/releases/download/' obj.release_tag '/' download_file];
        end

        function up_to_date = is_up_to_date(obj)
            up_to_date = strcmp(obj.installed_version, obj.upstream_version);
        end

        function success = install(obj)
            destination = [obj.install_path '.zip'];
            webopts = weboptions('Timeout', Inf);
            % check install or update
            installed_apps = ea_regexpdir(ea_prefsdir, '^SlicerForLeadDBS.*', 0, 'd');
            if isempty(installed_apps)
                msg_start = 'Installing';
            else
                warning('off','all')
                cellfun(@(x) rmdir(x, 's'), installed_apps);
                warning('on','all')
                msg_start = 'Updating';
            end
            f = msgbox([msg_start ' custom Slicer for Lead-DBS. This might take a while'], [msg_start '...'], 'help');
            % download
            try
                websave(destination, obj.download_url, webopts);
            catch
                delete(f);
                msgbox("Download error. Check you connection.", "Error", "error");
                success = 0;
                return
            end
            % install
            unzip(destination, fileparts(obj.install_path));
            delete(destination);

            % the zip file is packed with a subfolder inside it at the
            % moment, so we need to rename it.
            [~, extracted_folder] = fileparts(obj.download_url);
            extracted_folder = fullfile(ea_prefsdir, extracted_folder);
            try
                movefile(extracted_folder, obj.install_path);
                fclose(fopen(fullfile(obj.install_path, ['v' obj.upstream_version]), 'w'));
                success = 1;
                delete(f);
                msgbox("Custom Slicer for Lead-DBS installed!", "Installed", "help");
            catch ME
                success = 0;
                delete(f);
                msgbox(["Failed to install custom Slicer for Lead-DBS!", ME.message], "Failed", "Error");
            end
        end

        function update(obj)
            obj.install;
        end

        function status = run(obj, command)
            if ~exist('command','var')
                command='';
            end
            status = system([ea_path_helper(obj.executable_path) ' ' command]);
        end

    end
end

