function [] = ea_command_line_run(varargin)
% This function is called when lead has multiple inputs. These inputs are
% the patient directories and leadfigure handles with their values. Then,
% the specified command is executed via ea_run without a UI.
%
% For example, to run Segment normalization:
% lead dbs -normalize_checkbox -normmethod 2 /path/to/patient/directory

switch varargin{1}
    case {'dbs', '-d', 'd'}
        app = lead_dbs;
        leadprod = 'dbs';
    case {'connectome', 'conn', '-c', 'c'}
        app = lead_connectome;
        leadprod = 'connectome';
    case {'import', '-i', 'i'}
        leadprod = 'import';
end

% Import patient images to BIDS dataset
if strcmp(leadprod, 'import')
    images = varargin(isfile(varargin));

    % Set parameters
    varargin(isfile(varargin)) = [];
    for i = 2:2:length(varargin)
        switch lower(varargin{i}(2:end))
            case {'dataset', 'bids'}
                opt.dataset = varargin{i+1};
            case {'patientid', 'patient', 'id'}
                opt.patientid = varargin{i+1};
            case {'electrodemodel', 'ele', 'lead', 'model'}
                opt.electrodemodel = varargin{i+1};
        end
    end

    % Set electrode model to empty in case not specified
    if ~isfield(opt, 'electrodemodel')
        opt.electrodemodel = '';
    end

    ea_import_patient(opt.dataset, opt.patientid, images, opt.electrodemodel)
    return;
end

app.leadfigure.Visible = 'off';
handles = guidata(app.leadfigure);

% reset all checkboxes to 0
handles_names = fieldnames(handles);
for i = 1:length(handles_names)
    if isa(handles.(handles_names{i}), 'appdesigner.appmigration.UIControlPropertiesConverter') ...
            && strcmp(handles.(handles_names{i}).Style,'checkbox')
        handles.(handles_names{i}).Value = 0;
    end
end

patdirs = {};

for i = 2:nargin
    if isfolder(varargin{i})
        patdirs{end+1} = varargin{i};
    elseif startsWith(varargin{i}, '-')
        opt = varargin{i}(2:end);
        if strcmp(opt,'process')
            % run basic lead dbs pipeline with defaults
            handles.coreg_checkbox.Value = defaultOption('coreg_checkbox');
            handles.coregctmethod.Value = defaultOption('coregctmethod');
            handles.coregmrmethod.Value = defaultOption('coregmrmethod');
            handles.normalize_checkbox.Value = defaultOption('normalize_checkbox');
            handles.normmethod.Value = defaultOption('normmethod');
            handles.scrf.Value = defaultOption('scrf');
            handles.scrfmask.Value = defaultOption('scrfmask');
            handles.doreconstruction.Value = defaultOption('doreconstruction');
            handles.reconmethod.Value = defaultOption('reconmethod');
            handles.electrode_model_popup.Value = defaultOption('electrode_model_popup');
            handles.overwriteapproved.Value = defaultOption('overwriteapproved');
        elseif isfield(handles, opt)
            handles.(opt).Value = defaultOption(opt); % Set default option first as const value
        else
            error(['Unrecognized field: ' opt])
        end
    else % When the option value is specified explicitly
        if ismember(opt, {'coregctmethod', 'coregmrmethod', 'normmethod', 'scrfmask', 'reconmethod', 'electrode_model_popup'})
             if ~isnan(str2double(varargin{i})) % Number specified
                 handles.(opt).Value = str2double(varargin{i});
             else % Text specified
                Value = find(strcmp(handles.(opt).String, varargin{i}));
                if ~isempty(Value) % Text is exactly the value in the popupmenu
                    handles.(opt).Value = Value;
                else % Otherwise, do fuzzy match for some options
                    switch opt
                        case 'coregmrmethod'
                            if contains(varargin{i}, 'ANTs', 'IgnoreCase', 1)
                                handles.(opt).Value = 1;
                            elseif contains(varargin{i}, 'SPM', 'IgnoreCase', 1)
                                handles.(opt).Value = 8;
                            end
                        case 'normmethod'
                            if contains(varargin{i}, 'ANTs', 'IgnoreCase', 1)
                                handles.(opt).Value = 1;
                            elseif contains(varargin{i}, 'Apply', 'IgnoreCase', 1)
                                handles.(opt).Value = 2;
                            elseif contains(varargin{i}, 'EasyReg', 'IgnoreCase', 1)
                                handles.(opt).Value = 3;
                            elseif contains(varargin{i}, 'DARTEL', 'IgnoreCase', 1)
                                handles.(opt).Value = 6;
                            elseif contains(varargin{i}, 'Segment', 'IgnoreCase', 1)
                                handles.(opt).Value = 7;
                            elseif contains(varargin{i}, 'SHOOT', 'IgnoreCase', 1)
                                handles.(opt).Value = 8;
                            elseif contains(varargin{i}, 'SynthMorph', 'IgnoreCase', 1)
                                handles.(opt).Value = 9;
                            end
                        case 'scrfmask'
                            if contains(varargin{i}, 'No', 'IgnoreCase', 1)
                                handles.(opt).Value = 1;
                            elseif contains(varargin{i}, 'Coarse', 'IgnoreCase', 1)
                                handles.(opt).Value = 2;
                            elseif contains(varargin{i}, 'Fine', 'IgnoreCase', 1)
                                handles.(opt).Value = 3;
                            end
                        case 'reconmethod'
                            if contains(varargin{i}, 'Refined', 'IgnoreCase', 1)
                                handles.(opt).Value = 1;
                            elseif contains(varargin{i}, {'TRAC', 'CORE'}, 'IgnoreCase', 1)
                                handles.(opt).Value = 2;
                            elseif contains(varargin{i}, 'PaCER', 'IgnoreCase', 1)
                                handles.(opt).Value = 3;
                            elseif contains(varargin{i}, {'Manual', 'SPM'}, 'IgnoreCase', 1)
                                handles.(opt).Value = 4;
                            elseif contains(varargin{i}, 'Slicer', 'IgnoreCase', 1)
                                handles.(opt).Value = 5;
                            end
                    end
                end
             end
        else % For the CheckBox
            handles.(opt).Value = str2double(varargin{i});
        end
    end
end

options = ea_handles2options(handles);
options.leadprod = leadprod;
options.cmdlineCall = 1;
setappdata(app.leadfigure, 'handles', handles);
options.leadfigure = app.leadfigure;

% Erase last filesep
patdirs = erase(patdirs, filesep + textBoundary('end'));
% All patient folder should be under the same dataset
datasetDir = regexp(patdirs{1}, ['.*(?=\' filesep 'derivatives\' filesep 'leaddbs\' filesep 'sub-)'], 'match', 'once');
subjId = erase(patdirs, fullfile(datasetDir, 'derivatives', 'leaddbs', 'sub-'));
bids = BIDSFetcher(datasetDir);

options.uipatdirs = patdirs;
setappdata(options.leadfigure, 'bids', bids);
setappdata(options.leadfigure, 'subjId', subjId);

% dont show pop up methods
umachine = load(ea_prefspath('mat'));
machine = umachine.machine;
machine.methods_show = 0;
save(ea_prefspath('mat'), 'machine');

% run
ea_run('run',options);

% restore previous pop up config
machine.methods_show = umachine.machine.methods_show;
save(ea_prefspath('mat'), 'machine');

% Delete app
delete(app);


% Get default option values
function value = defaultOption(opt)
switch opt
    case 'coreg_checkbox'
        value = 1; % Do coregistration
    case 'coregctmethod'
        value = 1; % ANTs
    case 'coregmrmethod'
        value = 8; % SPM; value = 1 for ANTs
    case 'normalize_checkbox'
        value = 1; % Do normalization
    case 'normmethod'
        value = 1; % ANTs; 3 for EasyReg, 6 for SPM DARTEL, 7 for SPM Segment, 8 for SPM SHOOT, 9 for SynthMorph
    case 'scrf'
        value = 1; % Do brainshift correction
    case 'scrfmask'
        value = 2; % Coarse mask; 1 for No mask, 3 for Coarse + Fine Mask
    case 'doreconstruction'
        value = 1; % Do reconstruction
    case 'reconmethod'
        value = 1; % Refined TRAC/CORE; 2 for TRAC/CORE, 3 for PaCER, 4 for Manual, 5 for Slicer
    case 'electrode_model_popup'
        value = 1; % Medtronic 3389
    case 'overwriteapproved'
        value = 0;
end