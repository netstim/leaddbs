function ea_synctrajectoryhandles(handles,obj)
% set handles to match obj

% set coordinates
if obj.hasPlanning
set(handles.targetX,'String',num2str(obj.target.target(1))); set(handles.targetY,'String',num2str(obj.target.target(2))); set(handles.targetZ,'String',num2str(obj.target.target(3)));
set(handles.entryX,'String',num2str(obj.target.entry(1))); set(handles.entryY,'String',num2str(obj.target.entry(2))); set(handles.entryZ,'String',num2str(obj.target.entry(3)));
end
% set space
set(handles.space,'Value',obj.planRelative(5));

% set relative to radio buttons
switch obj.planRelative(1)
    case 1
        set(handles.AC,'Value',1);
    case 2
        set(handles.MCP,'Value',1);
    case 3
        set(handles.PC,'Value',1);
end
switch obj.planRelative(2)
    case 1
        set(handles.right,'Value',1);
    case 2
        set(handles.left,'Value',1);
end
switch obj.planRelative(3)
    case 1
        set(handles.anterior,'Value',1);
    case 2
        set(handles.posterior,'Value',1);
end
switch obj.planRelative(4)
    case 1
        set(handles.ventral,'Value',1);
    case 2
        set(handles.dorsal,'Value',1);
end

% set color backgroundcolor
set(handles.color,'BackgroundColor',obj.color);

%% set showplanning.
set(handles.showPlanning,'Value',obj.showPlanning);
set(handles.showPlanning,'enable',ea_bool2onoff(obj.hasPlanning));
% subordinate enables
if ~get(handles.showPlanning,'Value') || ~ea_bool2onoff(get(handles.showPlanning,'enable'))
    onoff='off';
else
    onoff='on';
end
set(handles.targetX,'enable',onoff); set(handles.targetY,'enable',onoff); set(handles.targetZ,'enable',onoff);
set(handles.entryX,'enable',onoff); set(handles.entryY,'enable',onoff); set(handles.entryZ,'enable',onoff);
set(handles.space,'enable',onoff);
set(handles.color,'enable',onoff);
switch obj.planningAppearance
    case 'line'
        set(handles.planningappearance,'Value',1);
        set(handles.electrode_model_plan,'Enable','off');
        set(handles.electrode_relative_plan,'Enable','off');
    case 'electrode'
        set(handles.planningappearance,'Value',2);
        set(handles.electrode_model_plan,'Enable','on');
        set(handles.electrode_relative_plan,'Enable','on');
end

string_list = get(handles.electrode_model_plan,'String');
[~,whichentry]=ismember(obj.plan2elstruct_model,string_list);
set(handles.electrode_model_plan,'Value',whichentry);


options.elmodel=obj.plan2elstruct_model;
options=ea_resolve_elspec(options);
useside=obj.side;
if useside>2
    useside=1;
end
if obj.electrodeRelativeToPlan>length(options.elspec.etagenames{useside})
    obj.electrodeRelativeToPlan=length(options.elspec.etagenames{useside});
elseif obj.electrodeRelativeToPlan==0 && options.elspec.tipiscontact
    obj.electrodeRelativeToPlan=1;
end
entrystring=options.elspec.etagenames{useside};
if ~options.elspec.tipiscontact
   entrystring=[{'Center of tip'},entrystring]; 
end
set(handles.electrode_relative_plan,'String',entrystring);

set(handles.electrode_relative_plan,'Value',obj.electrodeRelativeToPlan+double(~options.elspec.tipiscontact));

if ~(get(handles.showPlanning,'Value')) || ~ea_bool2onoff(get(handles.showPlanning,'enable')) || get(handles.space,'Value')>1
    onoff='off';
else
    onoff='on';
end
set(handles.right,'enable',onoff);
set(handles.left,'enable',onoff);
set(handles.anterior,'enable',onoff);
set(handles.posterior,'enable',onoff);
set(handles.ventral,'enable',onoff);
set(handles.dorsal,'enable',onoff);
set(handles.AC,'enable',onoff);
set(handles.PC,'enable',onoff);
set(handles.MCP,'enable',onoff);

%% macro/DBS electrode
set(handles.showMacro,'Value',obj.showMacro);
set(handles.showMacro,'enable',ea_bool2onoff(obj.hasMacro));
set(handles.electrode_model_popup,'String',ea_resolve_elspec);
% subordinate enables
if ~(get(handles.showMacro,'Value')) || ~ea_bool2onoff(get(handles.showMacro,'enable'))
    onoff='off';
else
    onoff='on';
end
set(handles.electrode_model_popup,'enable',onoff);

%% micro/MER
set(handles.showMicro,'Value',obj.showMicro);
set(handles.relateMicro,'String',{'Base location on Macroelectrode','Base location on Microelectrode'});
switch obj.relateMicro
    case 'macro'
        set(handles.relateMicro,'Value',1);
    case 'planning'
        set(handles.relateMicro,'Value',2);
end
% subordinate enables
set(handles.relateMicro,'enable',ea_bool2onoff(get(handles.showMicro,'Value')));
