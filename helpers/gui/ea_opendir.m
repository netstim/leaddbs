function ea_opendir(path, options)
% For single directory, change current working directory and open it in OS
% file manager.
%
% For single file, change current working directory to its parent folder
% and open the folder in OS file manager with the file selected.
%
% For multiple directories, open them in OS file manager.
%
% for multiple files, open their parent folders in OS file manager with the
% files selected.
%
% Do not change current working directory in case 'chdir' is set to false.

arguments
    path    {mustBeText}
    options.chdir   {mustBeNumericOrLogical} = true
end


if ischar(path)    % Single directory/file
    if ~isdeployed && options.chdir
        if isfolder(path)
            cd(path);
        elseif isfile(path)
            cd(fileparts(path));
        end
    end
    path = {path};
end

cellfun(@opendir, path);


function opendir(path)

path = GetFullPath(path);

if isfolder(path)
    if ismac
        system(['open "', path, '"']);
    elseif isunix
        system(['xdg-open "', path, '"']);
    elseif ispc
        system(['explorer "', path, '"']);
    end
elseif isfile(path)
    if ismac
        system(['open -R "', path, '"']);
    elseif isunix
        try
            system(['gdbus call --session --dest org.freedesktop.FileManager1 --object-path /org/freedesktop/FileManager1 --method org.freedesktop.FileManager1.ShowItems "[''file://', path, ''']" ""']);
        catch
            ea_cprintf('CmdWinWarnings', 'Failed to open directory:\n%s\n', fileparts(path));
        end
    elseif ispc
        system(['explorer /select,"', path, '"']);
    end
end
