function dirlist = ea_regexpdir(rootdir, expstr, recursive)
% wrapper for regexpdir (need to clear the persistent variable)

if ~exist('recursive','var')
    recursive = true;
end

% Fix to use '^STR' pattern recursively
if recursive && strcmp(expstr(1),'^')
    expstr = ['(^|.+[/\\])', expstr(2:end)];
end

clear regexpdir
dirlist = regexpdir(rootdir, expstr, recursive);
for i=1:length(dirlist)
    if contains(dirlist(i),'_Crop_1.nii');
        loc = strfind(dirlist{i},'_Crop_1.nii');
        dirlist{i} = strcat(dirlist{i}(1:loc-1),'.nii');
    end
end
    
