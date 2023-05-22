classdef ea_slicer_for_lead

    properties
        name
        install_path
        executable_path
    end
    
    methods
        function obj = ea_slicer_for_lead()
            if ismac
                os = 'macosx';
                date = '2022-11-03';
            elseif ispc
                os = 'win';
                date = '2022-11-03';
            else
                os = 'linux';
                date = '2023-01-20';                
            end
            obj.name = ['SlicerForLeadDBS-0.1.0-' date '-' os '-amd64'];
            obj.install_path = fullfile(ea_getearoot, 'ext_libs', obj.name);
            if ismac
                obj.executable_path = fullfile(obj.install_path, 'SlicerForLeadDBS.app', 'Contents', 'MacOS', 'SlicerForLeadDBS');
            elseif ispc
                obj.executable_path = fullfile(obj.install_path, 'SlicerForLeadDBS.exe');
            else
                obj.executable_path = fullfile(obj.install_path, 'SlicerForLeadDBS');
            end
        end


        function out = is_installed_and_up_to_date(obj)
            d = dir(fullfile(ea_getearoot, 'ext_libs', 'SlicerForLeadDBS*'));
            if isempty(d)
                out = 0;
                return
            end
            installed_version = regexp(d(1).name, '\d\d\d\d-\d\d-\d\d', 'match', 'once');
            latest_version = regexp(obj.name, '\d\d\d\d-\d\d-\d\d', 'match', 'once');
            out = strcmp(installed_version, latest_version);
        end
        
        function success = install(obj)
            destination = [obj.install_path '.zip'];
            downloadurl = ['https://github.com/netstim/SlicerForLeadDBS/releases/download/v221103-beta/' obj.name '.zip'];
            webopts = weboptions('Timeout', Inf);
            % check install or update
            installed_apps = dir(fullfile(ea_getearoot, 'ext_libs', 'SlicerForLeadDBS*'));
            if isempty(installed_apps)
                msg_start = 'Installing';
            else
                warning('off','all')
                cellfun(@(a,b) rmdir(fullfile(a,b),'s'), {installed_apps.folder}', {installed_apps.name}');
                warning('on','all')
                msg_start = 'Updating';
            end
            f = msgbox([msg_start ' custom Slicer for Lead-DBS. This might take a while'], [msg_start '...'], 'help');
            % download
            try
                websave(destination, downloadurl, webopts);
            catch
                delete(f);
                msgbox("Download error. Check you connection.", "Error", "error");
                success = 0;
                return
            end
            % install
            unzip(destination, fileparts(obj.install_path));
            delete(destination)
            success = 1;
            delete(f)
            msgbox("Custom Slicer for Lead-DBS installed!", "Installed", "help");
        end
        
        function [] = run(obj, command)
            if ~exist('command','var')
                command='';
            end
            system([obj.executable_path ' ' command]);
        end
        
    end
end

