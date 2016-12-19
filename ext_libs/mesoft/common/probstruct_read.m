function [res, errStr, fName]= probstruct_read(fNameIn, mFlag)
%
%function [res, errStr, fName]= probstruct_read(fNameIn, mFlag)
%
%   Opens one or a set of probstruct from the filesystem
%
%    mFlag == 'Single' (default)
%       fNameIn can be a filename or a path or empty (for file selection dialog). 
%       res is the opend probStrcut (converted to V2 if necessary) or is
%               empty (on error)
%       fName is a string containing the filename
%
%    mFlag == 'Multi' 
%       fNameIn can be a cell of filenames or a path or empty (for multi file selection dialog). 
%       res a cellarray containing the opended probStrcut (converted to V2
%           if necessary) or is empty (on error)
%
%
% Bjoern W. Kreher
% 06/08
%
% WinXP, UNIX

res= []; errStr= ''; fName= '';

if nargin < 2
    mFlag= 'Single';
end

if (nargin == 0)
    fNameIn = [];
end

%% jumb to defined directory (if necessary)
if ischar(fNameIn) && (exist(fNameIn, 'dir') == 7)
    curDir= pwd;
    cd(fNameIn);
    fNameIn = [];
else
    curDir= [];
end


% open file
[res, errStr, fName]= local_probstruct_read(fNameIn, mFlag);


%% jumb back to current directory (if necessary)
if ~isempty(curDir)
   cd(curDir);
end



%%%%%%%%5
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, fName]= local_probstruct_read(fNameIn, mFlag)
res= [];  errStr= ''; fName= [];
% file name selection
if strcmp(mFlag, 'Single')
    if isempty(fNameIn) || ~ischar(fNameIn)
        [fileStr,dirStr]=uigetfile('*.mat','Load a probMaps');
        if fileStr == 0
            errStr= sprintf('%s(warning): aborted by user', mfilename);
            return
        else
            fNameIn= fullfile(dirStr, fileStr);
        end
    end

    % check all fnames
    if ~ischar(fNameIn)
        errStr= sprintf('%s(error): file have to contain a strings', mfilename);
        return
    end
    
    [res, errStr]= local_openProbStruct(fNameIn);
elseif strcmp(mFlag, 'Multi')
    if isempty(fNameIn) || ~iscell(fNameIn)
        [fileStr,dirStr]=uigetfile('*.mat','Load one or more probMaps', 'MultiSelect', 'on');
        if ~iscell(fileStr) %one or none file selected
            if fileStr == 0
                errStr= sprintf('%s(warning): aborted by user', mfilename);
                return
            else
                fNameIn{1}= fullfile(dirStr, fileStr);
            end
        else
            fNameIn= cell(length(fileStr), 1);
            for i= 1:length(fileStr)
                fNameIn{i}= fullfile(dirStr, fileStr{i});
            end
        end
    end
    % check all fnames
    for i= 1:length(fNameIn)
        if ~ischar(fNameIn{i})
            errStr= sprintf('%s(error): cell array have to contain only strings', mfilename);
            return
        end
    end
    
    % open data
    res= cell(length(fNameIn), 1);
    for i= 1:length(fNameIn)
        [res{i}, errStr]= local_openProbStruct(fNameIn{i});        
        if isempty(res{i})
            res= [];
            return;
        end
    end
else
    errStr= sprintf('%s(error): invalid paramaters', mfilename);
    return
end

fName= fNameIn;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_openProbStruct(fName)
res= []; errStr= '';

if exist(fName,'file') ~= 2
    errStr= sprintf('%s::local_openProbStruct(error): File ''%s'' does not exist', ...
        mfilename, fName);
    return
end

mr= mrstruct_read(fName);
[ok, errStr, verStr]= probstruct_istype(mr);
if strcmp(verStr, 'V1')
    mr= probstruct_init(mr);
    if isempty(mr)
        errStr= sprintf('%s::local_openProbStruct(error): internal error', mfilename);
        return
    end
    [ok, errStr, verStr]= probstruct_istype(mr);
end

if ~ok
    errStr= sprintf('%s::local_openProbStruct(error): File ''%s'' is not a valid probStruct', mfilename, fName);
    res= [];
    return
end


res= mr;
