function saveSession(confirm)
% saveSession([confirm])
%
% confirm: ask before overwriting existing mrSession file
%          default: 0
%
% djh 5/2005

mrGlobals;

if ieNotDefined('confirm')
    confirm = 0;
end

pathStr = fullfile(MLR.homeDir,'mrSession.mat');

if confirm
    if exist(pathStr,'file')
        questionString = 'mrSession.mat already exists. Do you want to overwrite it?';
        buttonName = questdlg(questionString, 'Warning', 'Yes', 'No', 'No');
        pause(.1);  % Prevent hanging
        if strcmp(buttonName, 'No')
            return
        end
    end
else
    if strcmp(mrGetPref('verbose'),'Yes')
        disp('Warning: overwriting mrSession.mat');
    end
end

session = MLR.session;
groups = MLR.groups;
save(pathStr,'session','groups');
