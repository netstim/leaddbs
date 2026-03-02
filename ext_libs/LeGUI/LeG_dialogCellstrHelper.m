function outStr = LeG_dialogCellstrHelper (inputStr)
%   Copyright 2010-2011 The MathWorks, Inc.

% Helper used by MSGBOX, ERRORDLG, WARNDLG, QUESTDLG to parse the input
% string vector, matrix or cell array or strings.
% This works similar to the CELLSTR function but does not use deblank, like
% cellstr, to eliminate any trailing white spaces.

% Validate input string type. 
validateattributes(inputStr, {'char','cell'}, {'2d'},mfilename);

% Convert to cell array of strings without eliminating any user input. 
if ~iscell(inputStr)
    inputCell = {};
    for siz = 1:size(inputStr,1)
        inputCell{siz} =inputStr(siz,:); %#ok<AGROW>
    end
    outStr = inputCell;
else
    outStr = inputStr;
end
end
