function mrInitRet
% mrInitRet
%
% Creates mrSESSION data structure and Readme.txt file
%
% DJH, 11/2004

global MLR

% If mrSESSION.mat exists, prompt user to overwrite.
if exist('mrSESSION.mat','file')
    questionString = 'mrSESSION already exists. Do you want to continue, which will start from scratch?';
    buttonName = questdlg(questionString, 'Warning', 'Yes', 'No', 'No');
    pause(.1);  % Prevent hanging
    if strcmp(buttonName, 'No')
        return
    end
else
    MLR.homeDir = pwd;
end

% Initialize session
session = editSessionGUI;

% Initialize groups
groups = editGroupGUI('groupname','Raw');

% Save mrSessino
save mrSession session groups

% Create Readme.txt file
createReadme(session, groups);

clear all;
