function ea_importcorticalels(options)

% Import Cortical Electrodes

% Need to define input options

%% Left
start_dir = options.fsdir;
Side.Left = questdlg('Do you have a Left Sided Strip?');
if strcmp(Side.Left,'Yes')
    [CortElecFileName,CortElecFilePath,CortElecFileExt] = ea_uigetfile(start_dir,['Choose Cortical Electorde Coordinates File:  Left Side']);
    CortElecFile = fullfile(CortElecFilePath{1},[CortElecFileName{1},CortElecFileExt{1}]);
    CortElecStruct = load(CortElecFile);
    if isstruct(CortElecStruct(1))
        fields = fieldnames(CortElecStruct(1));
        for f = 1:length(fields);
            if iscell(eval(['CortElecStruct(1).' fields{f}])) && size(eval(['CortElecStruct(1).' fields{f}]),1)>size(eval(['CortElecStruct(1).' fields{f}]),2)
                CortElecs.Left = cell2mat(eval(['CortElecStruct(1).' fields{f}]));
            elseif iscell(eval(['CortElecStruct(1).' fields{f}])) && size(eval(['CortElecStruct(1).' fields{f}]),2)>size(eval(['CortElecStruct(1).' fields{f}]),1)
                CortElecs.Left = cell2mat(eval(['CortElecStruct(1).' fields{f}])');
            end
            % Check dim
            dims = size(CortElecs.Left);
            if isempty(dims(dims==3))
                ea_error('Error while importing cortical electrodes')
            end
        end
    end
elseif strcmp(Side.Left,'Cancel')
    return
end

%% Right
if exist('CortElecFilePath','var')
    start_dir = CortElecFilePath{1};
else
    start_dir = options.fsdir;
end

Side.Right = questdlg('Do you have a Right Sided Strip?');
if strcmp(Side.Right,'Yes')
    [CortElecFileName,CortElecFilePath,CortElecFileExt] = ea_uigetfile(start_dir,['Choose Cortical Electorde Coordinates File: Right Side']);
    CortElecFile = fullfile(CortElecFilePath{1},[CortElecFileName{1},CortElecFileExt{1}]);
    CortElecStruct = load(CortElecFile);
    if isstruct(CortElecStruct(1))
        fields = fieldnames(CortElecStruct(1));
        for f = 1:length(fields);
            if iscell(eval(['CortElecStruct(1).' fields{f}])) && size(eval(['CortElecStruct(1).' fields{f}]),1)>size(eval(['CortElecStruct(1).' fields{f}]),2)
                CortElecs.Right = cell2mat(eval(['CortElecStruct(1).' fields{f}]));
            elseif iscell(eval(['CortElecStruct(1).' fields{f}])) && size(eval(['CortElecStruct(1).' fields{f}]),2)>size(eval(['CortElecStruct(1).' fields{f}]),1)
                CortElecs.Right = cell2mat(eval(['CortElecStruct(1).' fields{f}])');
            end
            % Check dim
            dims = size(CortElecs.Left);
            if isempty(dims(dims==3))
                ea_error('Error while importing cortical electrodes')
            end
        end
    end
elseif strcmp(Side.Left,'Cancel')
    return
end

%% Other
% Side.Other = questdlg('Do you have any other Strips?');
% if strcmp(Side.Other,'Yes')
%     Side.Other_str = inputdlg('Please, Specify');
% end

% if strcmp(Side.Other,'Yes')
%     disp('Import Other Electrodes Feature unavailable.')
% end

%% Save To PatientDirectory/cortex
if exist('CortElecs','var')
disp(['Saving to ' options.uipatdirs,'/cortex/CortElecs.mat'])
    save([options.uipatdirs,'/cortex/CortElecs.mat'],'CortElecs')
end
