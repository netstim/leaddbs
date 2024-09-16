function [] = ea_command_line_run(varargin)
% This function is called when lead has multiple inputs. These inputs are
% the patient directories and leadfigure handles with their values. Then,
% the specified command is executed via ea_run without a UI.
%
% For example, to run Segment normalization:
% lead dbs -normalize_checkbox -normmethod 2 /path/to/patient/directory

switch varargin{1}
    case {'dbs', '-d', 'd'}
        h = lead_dbs;
        leadprod = 'dbs';
    case {'connectome', 'conn', '-c', 'c'}
        h = lead_connectome;
        leadprod = 'connectome';
end

h.leadfigure.Visible = 'off';
drawnow;
handles = guidata(h.leadfigure);

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
            handles.coreg_checkbox.Value = 1;
            handles.coregctmethod.Value = 1;
            handles.coregmrmethod.Value = 1;
            handles.normalize_checkbox.Value = 1;
            handles.normmethod.Value = 1;
            handles.scrf.Value = 1;
            handles.scrfmask.Value = 2;
            handles.doreconstruction.Value = 1;
            handles.reconmethod.Value = 3;
        elseif isfield(handles, opt)
            set(handles.(opt), 'Value', 1); % set default option to 1   
        else
            error(['Unrecognized field: ' opt])
        end
    else
        set(handles.(opt), 'Value', str2double(varargin{i}));
    end
end

options = ea_handles2options(handles);
options.leadprod = leadprod;
setappdata(h.leadfigure, 'handles', handles);
options.leadfigure = h.leadfigure;

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

h.leadfigure.Visible = 'on';
drawnow; % make gui visible so that close request function is executed
close(h)
