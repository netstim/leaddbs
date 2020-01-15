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

set(h, 'Visible', 'off'); drawnow;
handles = guihandles(h);

dirs = {};

for i = 2:nargin
    switch varargin{i}(1) % first character of str
        case filesep % directory
            patient_dir = varargin{i};
            if isfolder(patient_dir)
                dirs{end+1} = patient_dir;
            else
                error([patient_dir ' is not a folder'])
            end
        case '-' % option
            opt = varargin{i}(2:end);
            if isfield(handles, opt)
                set(handles.(opt),'Value',1); % set default option to 1
            elseif strcmp(opt,'process')
                % run basic lead dbs pipeline with defaults
                handles.coreg_checkbox.Value = 1;
                handles.normalize_checkbox.Value = 1;
                handles.scrf.Value = 1;
                handles.doreconstruction_checkbox.Value = 1;
            else
                error(['Unrecognized field: ' opt])
            end
        otherwise
            set(handles.(opt),'Value',str2double(varargin{i}));
    end
end

options = ea_handles2options(handles);
options.uipatdirs = dirs;
options.leadprod = leadprod;

% dont show pop up methods
umachine = load([ea_gethome, '.ea_prefs.mat']); 
machine = umachine.machine;
machine.methods_show = 0;
save([ea_gethome, '.ea_prefs.mat'],'machine');

% run
ea_run('run',options);

% restore previous pop up config
machine.methods_show = umachine.machine.methods_show;
save([ea_gethome, '.ea_prefs.mat'],'machine');


set(h, 'Visible', 'on'); drawnow; % make gui visible so that close request function is executed
close(h)
