% Mac OSX "quaranteens" downloaded executable files and prevents them from
% beeing executed for IT security reseasons. The user has to explicit remove
% the quranteen once, which is normally done by an explicit "right-click" and
% "open" user action. 
%
% This function checks if any executable files in the given folderpath are
% quranteend and thus cannot be executed preventing software reling on them
% from funtioning properly. 
%
% If quranteened files are found, a warning is raised an a solution command
% to lift the file qurantine is suggested to the user.
%
% Andreas Husch, andreas.husch@uni.lu
% University of Luxembourg, LCSB, Imaging AI Group, 2023

function checkMacOSFileQuarantine(folderpath, filepattern)
    if(nargin < 2)
        filepattern = '*.mexmaci64';
    end

    [~, xattrs] = system(['xattr ' folderpath filesep filepattern]);

    % check for files that are in status 'com.apple.quarantine'
    xattrsCell = split(xattrs, newline);
    quarantinedStr = strfind(xattrsCell, 'com.apple.quarantine');
    qurantinedIdxs = cellfun(@(cell)(~isempty(cell)), quarantinedStr);

    if(any(qurantinedIdxs))
        warning(['Found executable ' filepattern ' files in folder ' folderpath ' that are under Mac OSX file quarantine and thus cannot be executed: ' newline newline ...
            xattrsCell{qurantinedIdxs} newline newline...
            'lift the qurantine by "right-click" and "open" each file once or by running terminal command' newline ...
            '"sudo xattr -r -d com.apple.quarantine ' folderpath filesep filepattern '"'] )
    end
end