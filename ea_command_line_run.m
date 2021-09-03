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

% reset all checkboxes to 0
handles_names = fieldnames(handles);
for i = 1:length(handles_names)
    if isprop(handles.(handles_names{i}), 'Style') & strcmp(handles.(handles_names{i}).Style,'checkbox')
        handles.(handles_names{i}).Value = 0;
    end
end

tryBIDS = false; % flag to try BIDS import when calling with -process
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
                handles.coregmrmethod.Value = 1;
                handles.normalize_checkbox.Value = 1;
                handles.normmethod.Value = 9;
                handles.scrf.Value = 1;
                handles.scrfmask.Value = 2;
                handles.doreconstruction_checkbox.Value = 1;
                handles.reconmethod.Value = 3;
                tryBIDS = true;
            else
                error(['Unrecognized field: ' opt])
            end
        otherwise
            set(handles.(opt),'Value',str2double(varargin{i}));
    end
end


if tryBIDS
    try 
        for i = 1:length(dirs) % BIDS subject dir
            spm_BIDS(fileparts(fileparts(fullfile(dirs{i},filesep))));
        end
        handles.dicomcheck.Value = 1;
        handles.dcm2niiselect.Value = 4;
    end
    if length(dirs) == 1 % root BIDS directory
        try
            spm_BIDS(dirs{1});
            D = dir(fullfile(dirs{1},'sub*'));
            dirs = cellstr(strcat(vertcat(D.folder),repmat(filesep,length(D),1),vertcat(D.name)));
            handles.dicomcheck.Value = 1;
            handles.dcm2niiselect.Value = 4;
        end
    end
end

options = ea_handles2options(handles);
options.uipatdirs = dirs;
options.leadprod = leadprod;
setappdata(h,'handles',handles);
options.leadfigure = h;

% dont show pop up methods
umachine = load([ea_gethome, '.ea_prefs.mat']); 
machine = umachine.machine;
machine.methods_show = 0;
save([ea_gethome, '.ea_prefs.mat'],'machine');

% run
ea_run('run',options);

if handles.dicomcheck.Value == 1 && tryBIDS % first run only imports bids folders. run again to process 
    handles.dicomcheck.Value = 0;
    options = ea_handles2options(handles);
    options.uipatdirs = strsplit(handles.patdir_choosebox.Tooltip);    
    options.leadprod = leadprod;
    ea_run('run',options);
end

% restore previous pop up config
machine.methods_show = umachine.machine.methods_show;
save([ea_gethome, '.ea_prefs.mat'],'machine');


set(h, 'Visible', 'on'); drawnow; % make gui visible so that close request function is executed
close(h)
