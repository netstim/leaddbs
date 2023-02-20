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
                date = '2022-11-03'
            elseif ispc
                os = 'win';
                date = '2022-11-03'
            else
                os = 'linux';
                date = '2023-01-20'
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
        
        function out = is_installed(obj)
            out = isfile(obj.executable_path);
        end
        
        function success = install(obj)
            destination = [obj.install_path '.zip'];
            downloadurl = ['https://github.com/netstim/SlicerForLeadDBS/releases/download/v221103-beta/' obj.name '.zip'];
            webopts = weboptions('Timeout', Inf);
            f = msgbox("Installing custom Slicer for Lead-DBS. This might take a while", "Installing...", "help");
            try
                websave(destination, downloadurl, webopts);
            catch
                delete(f);
                msgbox("Download error. Check you connection.", "Error", "error");
                success = 0;
                return
            end
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

