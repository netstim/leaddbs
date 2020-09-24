function ea_opendir(directory)
% For single directory, change current directory to it and open it in OS
% file manager. For multiple directories, open them in OS file manager.

if ischar(directory)    % Single directory
    if ~isdeployed
        cd(directory);
    end
    directory = {directory};
end

for d=1:length(directory)
    if ismac
        system(['open "', directory{d}, '"']);
    elseif isunix
        system(['xdg-open "', directory{d}, '"']);
    elseif ispc
        system(['explorer "', directory{d}, '"']);
    end
end
