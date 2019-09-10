function [] = ea_command_line_run(varargin)
% This function is called when lead has multiple inputs. These inputs are
% the patient directories and leadfigure handles with their values. Then,
% the specified command is executed via ea_run without a UI.
%
% For example, to run Segment normalization:
% lead dbs /path/to/patient/directory -normalize_checkbox -normmethod 2


switch varargin{1}
    case 'dbs'
        h = lead_dbs;
    case 'connectome'
        h = lead_connectome;
end

set(h, 'Visible', 'off'); pause(1);
handles = guihandles(h);

dirs = {};

for i = 2:nargin
    switch varargin{i}(1) % first character of str
        case filesep % directory
            dirs{end+1} = varargin{i};
        case '-' % option
            if i+1 <= nargin && varargin{i+1}(1) ~= '-' % check if the option is set with value
                set(handles.(varargin{i}(2:end)),'Value',str2double(varargin{i+1}));
                i = i+1;
            else
                set(handles.(varargin{i}(2:end)),'Value',1); % set default option to 1
            end     
    end
end

options = ea_handles2options(handles);
options.uipatdirs = dirs;
options.leadprod = varargin{1};

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


set(h, 'Visible', 'on'); pause(1); % make gui visible so that close request function is executed
close(h)
