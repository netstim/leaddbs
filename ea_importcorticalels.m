function CortElecs = ea_importcorticalels(options)

% Import Cortical Electrodes
% __________________________________________________________________________________
% Copyright (C) 2017 University of Pittsburgh (UPMC), Brain Modulation Lab
% Ari Kappel

% save info to struct
CortElecs.patientname = options.patientname;
CortElecs.ptdir = options.uipatdirs;
CortElecs.fsdir = options.fsdir;

%% Determine Side
start_dir = options.fsdir;
eSide = ea_popupquest('On which side are your electrodes located?',...
        'Right','Left','Both');

% Left or Right    
switch eSide
    case {'Left','Right'}
        [CortElecFileName,CortElecFilePath,CortElecFileExt] = ea_uigetfile(start_dir,['Choose Cortical Electorde Coordinates File:  ' eSide ' Side']);
        CortElecFile = fullfile(CortElecFilePath{1},[CortElecFileName{1},CortElecFileExt{1}]);
        [CortElecs.(eSide).pathname,file,ext] = fileparts(CortElecFile);
        CortElecs.(eSide).filename = [file,ext];
        CortElecs.(eSide).nativespace = ea_popupquest('In what space were your electrode locations determined?',...
            'Preop','Postop');
        
    CortElecStruct = load(CortElecFile);
    if isstruct(CortElecStruct(1))
        fields = fieldnames(CortElecStruct(1));
        for sMRI = 1:length(fields);
            if iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),1)>size(eval(['CortElecStruct(1).' fields{sMRI}]),2)
                CortElecs.(eSide).raw = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}]));
            elseif iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),2)>size(eval(['CortElecStruct(1).' fields{sMRI}]),1)
                CortElecs.(eSide).raw = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}])');
            end
            % Check dim
            dims = size(CortElecs.(eSide).raw);
            if size(CortElecs.(eSide).raw,2)~=3
                ea_error('Error while importing cortical electrodes')
            end
        end
    end
    % Reading in MRI parameters
    sMRI=MRIread(fullfile(options.fsdir,'mri/T1.nii'));
    
    coregfile = fdir(ptdir,'2postop');
    % ONLY SUPPORTS ANTS COREGISTRATION
    if ~isempty(coregfile) && length(coregfile)==1 && ~isempty(regexp(coregfile,'\wants','match')) 
        coregfile = cell2mat(coregfile);
        load([ptdir,filesep,coregfile])
        aff = ea_antsmat2mat(AffineTransform_float_3_3,fixed);
    end
        aff([4,8,12])=aff([13:15]); aff([13:15])=0;
    % Translating into the appropriate space
    for k=1:size(CortElecs.(eSide).raw,1)
        a=sMRI.vox2ras/sMRI.tkrvox2ras*[CortElecs.(eSide).raw(k,:) 1]';
        CortElecs.(eSide).native(k,:)=a(1:3)';
    end
    CortElecs.(eSide).mni = cs_convert(sMRI, 'scs', 'mni', CortElecs.Left.native);
    

% Both
case 'Both'
    cellSide = {'Right','Left'};
    for i = 1:2
        [CortElecFileName,CortElecFilePath,CortElecFileExt] = ea_uigetfile(start_dir,['Choose Cortical Electorde Coordinates File:  Left Side']);
        CortElecFile = fullfile(CortElecFilePath{1},[CortElecFileName{1},CortElecFileExt{1}]);
        [CortElecs.(cellSide{i}).pathname,file,ext] = fileparts(CortElecFile);
        CortElecs.(cellSide{i}).filename = [file,ext];
        CortElecs.(cellSide{i}).nativespace = ea_popupquest('In what space were your electrode locations determined?',...
            'Preop','Postop');
        
    CortElecStruct = load(CortElecFile);
    if isstruct(CortElecStruct(1))
        fields = fieldnames(CortElecStruct(1));
        for sMRI = 1:length(fields);
            if iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),1)>size(eval(['CortElecStruct(1).' fields{sMRI}]),2)
                CortElecs.(cellSide{i}).raw = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}]));
            elseif iscell(eval(['CortElecStruct(1).' fields{sMRI}])) && size(eval(['CortElecStruct(1).' fields{sMRI}]),2)>size(eval(['CortElecStruct(1).' fields{sMRI}]),1)
                CortElecs.(cellSide{i}).raw = cell2mat(eval(['CortElecStruct(1).' fields{sMRI}])');
            end
            % Check dim
            dims = size(CortElecs.(cellSide{i}).raw);
            if size(CortElecs.(cellSide{i}).raw,2)~=3
                ea_error('Error while importing cortical electrodes')
            end
        end
    end
    % Reading in MRI parameters
    sMRI=MRIread(fullfile(options.fsdir,'mri/T1.nii'));
    
    % Translating into the appropriate space
    for k=1:size(CortElecs.(cellSide{i}).raw,1)
        a=sMRI.vox2ras/sMRI.tkrvox2ras*[CortElecs.(cellSide{i}).raw(k,:) 1]';
        CortElecs.(cellSide{i}).native(k,:)=a(1:3)';
    end
    CortElecs.(cellSide{i}).mni = cs_convert(sMRI, 'scs', 'mni', CortElecs.Left.native);
    end 
end


function files = fdir(Path,exp)
% Returns contents of path as cell of filenames

if isempty(Path) || ~exist('Path','var')
    Path = pwd;
end
    
contents = dir(Path);
contents = contents(cellfun(@(x) isempty(regexp(x, '^\.', 'once')), {contents.name}));
files = {contents(~[contents.isdir]).name};

if nargin==2
    files = files(~cellfun(@(x) isempty(regexp(x, ['\w' exp],'match')), files));
end

