% Converts all .m and .mat files in the current directory to lowercase filenames
%		move_files_lowercase
%
%
%		Examples:
%		move_files_lowercase;
%
%		Harald Fischer
%		1/00
%
%		PC



function files_lowercase_rec(mPath)

moveCall= 'mv';
mPath
%%% handle m-files
dirStruct           = dir(mPath);
[noOfEntries dummy] = size(dirStruct);
directories= {};
for k=1:noOfEntries
    fileStr = dirStruct(k).name;
    if dirStruct(k).isdir==1 && strcmp(fileStr,'..')==0 && ...
            strcmp(fileStr,'.')==0
        directories{length(directories)+1}= lower(fileStr);
        cmdString = sprintf('%s %s/%s %s/%s', moveCall, ...
            mPath, fileStr, ...
            mPath, lower(fileStr));
        if strcmp(fileStr, lower(fileStr)) == 0
            unix(cmdString);
        end
    elseif (strcmp(fileStr,'..')==0 && strcmp(fileStr,'.')==0) && ...
            (strcmp(fileStr(length(fileStr)-1:length(fileStr)), '.M') || ...
            strcmp(fileStr(length(fileStr)-1:length(fileStr)), '.m') || ...
            strcmp(fileStr(length(fileStr)-3:length(fileStr)), '.FIG') || ...
            strcmp(fileStr(length(fileStr)-3:length(fileStr)), '.fig') || ...
            strcmp(fileStr(length(fileStr)-3:length(fileStr)), '.MAT') || ...
            strcmp(fileStr(length(fileStr)-3:length(fileStr)), '.mat'))
        cmdString = sprintf('%s %s/%s %s/%s', moveCall, ...
            mPath, fileStr, ...
            mPath, lower(fileStr));
        if strcmp(fileStr, lower(fileStr)) == 0
            unix(cmdString);
        end
    end
end
%%% End of: handle m-files

for k=1:length(directories)
    if strcmp(directories{k}, 'cvs') == 0 
        newPath= sprintf('%s/%s', mPath, directories{k});
        files_lowercase_rec(newPath);
    end
end

% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics