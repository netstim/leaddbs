function [nameStr, errStr]= get_unique_str(entriesCell, defaultStr, nameGenFlag, strucFlag, forceDialogFlag)
%function [nameStr, errStr]= get_unique_str(entriesCell, defaultStr, nameGenFlag, strucFlag, forceDialogFlag)
%
% 
% Bjoern W. Kreher
% 03/04
%
% UNIX


% Set default flags
if ~exist('defaultStr') || isempty(defaultStr)
    defaultStr= 'default';
end
if ~exist('nameGenFlag') || isempty(nameGenFlag)
    nameGenFlag= 1;
end
if ~exist('strucFlag') || isempty(strucFlag)
    strucFlag= 1;
end
if ~exist('forceDialogFlag') || isempty(forceDialogFlag)
    forceDialogFlag= 0;
end
nameStr= ''; errStr= ''; ok= 0;


if nameGenFlag
    tmpStr= local_genNewUniqueName(entriesCell, defaultStr, '%s_%03d');
else
    tmpStr= defaultStr;   
end
%nameStr= tmpStr(1:(length(tmpStr) + 1 - max(find(cumsum(abs([' ' tmpStr(end:-1:1)]-32)) == 0)))); % eleminate successing spaces
nameStr= tmpStr;

% Condition in case dialog are not forced
if ~forceDialogFlag
    [ok, errStr]= local_testString(entriesCell, nameStr, strucFlag);
    if ok
        return;
    else
        uiwait(warndlg(errStr, 'Invalid Name', 'modal'))
    end
end

% mail loop
while ~ok
    nameStrCell= inputdlg('Please enter a unique name', 'Name selection', 1, {nameStr});
    if isempty(nameStrCell)
        nameStr= ''; 
        errStr= 'Input was aborted by user!';
        return;
    else
        tmpStr= nameStrCell{1};
        %nameStr= tmpStr(1:(length(tmpStr) + 1 - max(find(cumsum(abs([' ' tmpStr(end:-1:1)]-32)) == 0)))); % eleminate successing spaces
        nameStr= tmpStr;
    end
    [ok, errStr]= local_testString(entriesCell, nameStr, strucFlag);
    if ~ok
        uiwait(warndlg(errStr, 'Invalid Name', 'modal'))
    end
end


%
%
%  local Funktions
%
%

% function test if 'nameStr' is unique in entriesCell. if not, errStr containes message
function [ok, errStr]= local_testString(entriesCell, nameStr, strucConstraintFlag)

if isempty(nameStr) || any(strcmp(entriesCell, nameStr))
    errStr=strcat('The name ''', nameStr, ''' is not unique! Pleas try an other name.');
    ok= 0;
elseif strucConstraintFlag && ~local_isfield_compatible(nameStr)
    errStr=strcat('The name ''', nameStr, ''' contains invalid characters. Pleas don''t use +, -, *, /,  , etc.');
    ok= 0;
else
    errStr= '';
    ok= 1;    
end


% function generates unique name by merging the 'defaultStr' and a counter
function nameStr= local_genNewUniqueName(entriesCell, defaultStr, patterStr)

if ~exist('patterStr') || isempty(patterStr)
    patterStr= '%s_%03d';    
end
ok= 0; i= 1;

while ~ok
    nameStr= sprintf(patterStr, defaultStr, i);
    ok= ~any(strcmp(entriesCell, nameStr));
    i= i + 1;
end


%test funktion ob 'nameStr' als feldname verwendet werden kann
function ok= local_isfield_compatible(nameStr)

ok= 1;
if any(nameStr == ' ')
    ok= 0;
    return
end
try cell2struct({[]}, nameStr);
catch, ok= 0;
end