function CortElecs = ea_importcorticalels(options)

% Import Cortical Electrodes
% Need to define input options

% future compatibility
patientname = options.patientname;
ptdir = options.ptdir;
fsdir = options.fsdir;

% save info to struct
CortElecs.patientname = patientname;
CortElecs.ptdir = ptdir;
CortElecs.fsdir = fsdir;
%% Left
start_dir = options.fsdir;
Side.Left = questdlg('Do you have a Left Sided Strip?');
if strcmp(Side.Left,'Yes')
    [CortElecFileName,CortElecFilePath,CortElecFileExt] = ea_uigetfile(start_dir,['Choose Cortical Electorde Coordinates File:  Left Side']);
    CortElecFile = fullfile(CortElecFilePath{1},[CortElecFileName{1},CortElecFileExt{1}]);
    CortElecs.LeftCortElecFile = CortElecFile;
    CortElecStruct = load(CortElecFile);
    if isstruct(CortElecStruct(1))
        fields = fieldnames(CortElecStruct(1));
        for sMRI = 1:length(fields);
            if iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),1)>size(eval(['CortElecStruct(1).' fields{sMRI}]),2)
                CortElecs.Left_raw = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}]));
            elseif iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),2)>size(eval(['CortElecStruct(1).' fields{sMRI}]),1)
                CortElecs.Left_raw = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}])');
            end
            % Check dim
            dims = size(CortElecs.Left_raw);
            if size(CortElecs.Left_raw,2)~=3
                ea_error('Error while importing cortical electrodes')
            end
        end
    end
    % Reading in MRI parameters
    sMRI=MRIread(fullfile(options.fsdir,'mri/T1.nii'));
    
    % Translating into the appropriate space
    for k=1:size(CortElecs.Left_raw,1)
        a=sMRI.vox2ras/sMRI.tkrvox2ras*[CortElecs.Left_raw(k,:) 1]';
        CortElecs.Left_native(k,:)=a(1:3)';
    end
    CortElecs.Left_mni = cs_convert(sMRI, 'scs', 'mni', CortElecs.Left_native)
    
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
    CortElecs.RightCortElecFile = CortElecFile;
    CortElecStruct = load(CortElecFile);
    if isstruct(CortElecStruct(1))
        fields = fieldnames(CortElecStruct(1));
        for sMRI = 1:length(fields);
            if iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),1)>size(eval(['CortElecStruct(1).' fields{sMRI}]),2)
                CortElecs.Right = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}]));
            elseif iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),2)>size(eval(['CortElecStruct(1).' fields{sMRI}]),1)
                CortElecs.Right = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}])');
            end
            % Check dim
            dims = size(CortElecs.Left);
            if isempty(dims(dims==3))
                ea_error('Error while importing cortical electrodes')
            end
        end
    end
elseif strcmp(Side.Right,'Cancel')
    return
end
