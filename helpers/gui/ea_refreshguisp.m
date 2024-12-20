function ea_refreshguisp(varargin)

handles=varargin{1};
options=varargin{2};

elstruct=getappdata(handles.stimfig,'elstruct');
S=getappdata(handles.stimfig,'S');
groupmode=getappdata(handles.stimfig,'groupmode');
actpt=getappdata(handles.stimfig,'actpt');

if isempty(actpt) || length(actpt)>1
    actpt=1;
end

if groupmode
    grouploaded=getappdata(handles.stimfig,'grouploaded');
    if isempty(grouploaded) % this is done only once and gets the selection info from lead_group initially (which patient shown).
        lgfig=getappdata(handles.stimfig,'resultfig');
        M=getappdata(lgfig,'M');
        actpt=M.ui.listselect;

        if length(actpt)>1 % more than one entry selected
            actpt=1;
        end

        % Ensure active patient is non empty
        % This can happen if you delete a patient, then add a new one, without clicking on the patient window
        if isempty(actpt)
            actpt=1;
        end
        setappdata(handles.stimfig,'actpt',actpt);
        % set grouploaded true is being done below.
    end

    elstruct=getappdata(handles.stimfig,'elstruct');

    set(handles.headertxt,'String',['Patient (',num2str(actpt),'/', num2str(length(elstruct)),'): ',elstruct(actpt).name]);

    gSv=getappdata(handles.stimfig,'gSv');
    if isfield(gSv,'vatmodel')
        if isempty(gSv.vatmodel)
            nms=get(handles.modelselect,'String');
            try
                gSv.vatmodel=nms{get(handles.modelselect,'Value')};
            catch
                keyboard
            end
            setappdata(handles.stimfig,'gSv',gSv);
        else
            [~,ind]=ismember(gSv.vatmodel,get(handles.modelselect,'String'));
            set(handles.modelselect,'Value',ind);
        end
    else
        nms=get(handles.modelselect,'String');
        try
            gSv.vatmodel=nms{get(handles.modelselect,'Value')};
        catch
            keyboard
        end
        setappdata(handles.stimfig,'gSv',gSv);
    end

    % load gS - updated with each refresh:
    gS=getappdata(handles.stimfig,'gS');
    if isempty(grouploaded)
        if ~isempty(gS)
            % determine stimlabel from priorly set gS:
            for sub=1:length(gS)
                stimlabel=['gs_',M.guid];
                if ~isempty(stimlabel)
                    break
                end
            end
            setappdata(handles.stimfig,'stimlabel',stimlabel);

            % if gS is defined but group has just now been loaded
            try
                if ~isempty(gS(actpt).Rs1) % current patient is defined -> set S to gS of this patient.
                    S = gS(actpt);
                end
            end
        end

        % now tell everyone that the figure has been opened for a while already:
        setappdata(handles.stimfig,'grouploaded',1);
    end
end

stimlabel=getappdata(handles.stimfig,'stimlabel');

if isempty(S)
    S = ea_initializeS(stimlabel,options,handles);
    if ~isempty(S)
        try
            S.model = gSv.vatmodel;
        end
        setappdata(handles.stimfig, 'stimlabel', S.label);
    else
        setappdata(handles.stimfig, 'Status', 'Cancelled');
        return;
    end
else
    if isempty(S.Rs1)
        S = ea_initializeS(stimlabel,options,handles);
        if ~isempty(S)
            try
                S.model = gSv.vatmodel;
            end
            setappdata(handles.stimfig, 'stimlabel', S.label);
        else
            setappdata(handles.stimfig, 'Status', 'Cancelled');
            return;
        end
    end
end

if isfield(S, 'model')
    [~,ix]=ismember(S.model,get(handles.modelselect,'String'));
    if ~ix
        ea_error('The model of the selected stimulation is not available.');
    else
        set(handles.modelselect,'Value',ix);
        model = handles.modelselect.String{ix};
    end
else
    set(handles.modelselect,'Value',1);
    model = handles.modelselect.String{1};
end

vatsettings.butenko_calcPAM = ea_getprefs('vatsettings.butenko_calcPAM');
% if contains(model, 'OSS-DBS') && vatsettings.butenko_calcPAM
%     handles.Rs2am.Value = 0;
%     handles.Rs3am.Value = 0;
%     handles.Rs4am.Value = 0;
%     handles.Ls2am.Value = 0;
%     handles.Ls3am.Value = 0;
%     handles.Ls4am.Value = 0;
%     S.active = [1 1];
% end

Ractive=S.active(1);
Lactive=S.active(2);

if nargin==3
    if ischar(varargin{3})
        if endsWith(varargin{3}, 'R')
            side = 1;
        elseif endsWith(varargin{3}, 'L')
            side = 2;
        end
        sidestr = varargin{3}(end);

        if startsWith(varargin{3}, 'C', 'IgnoreCase', true) % CaseR or CaseL
            S=ea_redistribute_voltage(S,varargin{3});
            if S.([sidestr,'s',num2str(S.active(side))]).case.pol==1
                S.([sidestr,'s',num2str(S.active(side))]).case.pol=2;
                S=ea_redistribute_voltage(S,varargin{3});
                for k=1:S.numContacts
                    if S.([sidestr,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==2
                        S.([sidestr,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol=0;
                        S=ea_redistribute_voltage(S,['k',num2str(k),sidestr]);
                    end
                end
            end
        elseif startsWith(varargin{3}, 'k') % k#R or K#L
            S=ea_redistribute_voltage(S,varargin{3});
            if S.([sidestr,'s',num2str(S.active(side))]).(varargin{3}(1:end-1)).pol==2 && S.([sidestr,'s',num2str(S.active(side))]).case.pol==2
                S.([sidestr,'s',num2str(S.active(side))]).case.pol=0;
                S=ea_redistribute_voltage(S,['Case',sidestr]);
            end
        end
    else
        S=ea_redistribute_voltage(S,varargin{3});
    end
end

if contains(model, 'OSS-DBS')
    vatsettings.butenko_pulseWidth = ea_getprefs('vatsettings.butenko_pulseWidth');
    for source = 1:4
        if ~isfield(S.(['Rs', num2str(source)]), 'pulseWidth')
            S.(['Rs', num2str(source)]).pulseWidth=vatsettings.butenko_pulseWidth;
        end
        if ~isfield(S.(['Ls', num2str(source)]), 'pulseWidth')
            S.(['Ls', num2str(source)]).pulseWidth=vatsettings.butenko_pulseWidth;
        end
    end
else
    for source = 1:4
        if isfield(S.(['Rs', num2str(source)]), 'pulseWidth')
            S.(['Rs', num2str(source)]) = rmfield(S.(['Rs', num2str(source)]), 'pulseWidth');
        end
        if isfield(S.(['Ls', num2str(source)]), 'pulseWidth')
            S.(['Ls', num2str(source)]) = rmfield(S.(['Ls', num2str(source)]), 'pulseWidth');
        end
    end
end

setappdata(handles.stimfig,'S',S);

% set stim amplitudes
for source=1:4
    S.amplitude{1}(source)=S.(['Rs',num2str(source)]).amp;

    set(handles.(['Rs',num2str(source),'am']),'String',num2str(S.amplitude{1}(source)));
    set(handles.(['Rs',num2str(source),'va']),'Value',S.(['Rs',num2str(source)]).va);

    %if S.(['Rs',num2str(source)]).amp % check if a valid +/- combination is active, if not set defaults.
    anycontactpositive=0; anycontactnegative=0;
    for k=1:S.numContacts
        if S.(['Rs',num2str(source)]).(['k',num2str(k)]).pol==1
            anycontactnegative=1;
        elseif S.(['Rs',num2str(source)]).(['k',num2str(k)]).pol==2
            anycontactpositive=2;
        end
    end

    if ~anycontactnegative
        S.(['Rs',num2str(source)]).k2.pol=1;
        S.(['Rs',num2str(source)]).k2.perc=100;
    end

    if ~anycontactpositive
        S.(['Rs',num2str(source)]).case.pol=2;
        S.(['Rs',num2str(source)]).case.perc=100;
    end
    %end
end

for source=1:4
    S.amplitude{2}(source)=S.(['Ls',num2str(source)]).amp;

    set(handles.(['Ls',num2str(source),'am']),'String',num2str(S.amplitude{2}(source)));
    set(handles.(['Ls',num2str(source),'va']),'Value',S.(['Ls',num2str(source)]).va);

    % if S.(['Ls',num2str(source)]).amp % check if a valid +/- combination is active, if not set defaults.
    anycontactpositive=0; anycontactnegative=0;
    for k=1:S.numContacts
        if S.(['Ls',num2str(source)]).(['k',num2str(k)]).pol==1
            anycontactnegative=1;
        elseif S.(['Ls',num2str(source)]).(['k',num2str(k)]).pol==2
            anycontactpositive=1;
        end
    end

    if ~anycontactnegative
        S.(['Ls',num2str(source)]).k2.pol=1;
        S.(['Ls',num2str(source)]).k2.perc=100;
    end
    if ~anycontactpositive
        S.(['Ls',num2str(source)]).case.pol=2;
        S.(['Ls',num2str(source)]).case.perc=100;
    end
    % end
end

%% model to handles: all GUI elements.
source=Ractive;
for k=1:S.numContacts
    val=S.(['Rs',num2str(source)]).(['k',num2str(k)]).perc;
    set(handles.(['k',num2str(k),'uR']),'String',num2str(val));

    val=S.(['Rs',num2str(source)]).(['k',num2str(k)]).imp;
    set(handles.(['k',num2str(k),'imR']),'String',num2str(val));
end

% set case
set(handles.CuR,'String',num2str(S.(['Rs',num2str(source)]).case.perc));

source=Lactive;
for k=1:S.numContacts
    val=S.(['Ls',num2str(source)]).(['k',num2str(k)]).perc;
    set(handles.(['k',num2str(k),'uL']),'String',num2str(val));

    val=S.(['Ls',num2str(source)]).(['k',num2str(k)]).imp;
    set(handles.(['k',num2str(k),'imL']),'String',num2str(val));
end

% set case
set(handles.CuL,'String',num2str(S.(['Ls',num2str(source)]).case.perc));

if contains(model, 'OSS-DBS')
    handles.pulseWidthTextbox_R.String = num2str(S.(['Rs',num2str(S.active(1))]).pulseWidth);
    handles.pulseWidthTextbox_L.String = num2str(S.(['Ls',num2str(S.active(2))]).pulseWidth);
end

%% model to handles: polarity toggles
for k=1:S.numContacts
    if S.(['Rs',num2str(Ractive)]).(['k',num2str(k)]).pol==0 % off
        icon=fullfile(options.earoot, 'icons', ['empty',num2str(Ractive),'.png']);
    elseif S.(['Rs',num2str(Ractive)]).(['k',num2str(k)]).pol==1 % negative S1
        icon=fullfile(options.earoot, 'icons', ['minus',num2str(Ractive),'.png']);
    elseif S.(['Rs',num2str(Ractive)]).(['k',num2str(k)]).pol==2 % positive S1
        icon=fullfile(options.earoot, 'icons', ['plus',num2str(Ractive),'.png']);
    end
    handles.(['k',num2str(k),'polR']).ImageSource = icon;
    handles.(['k',num2str(k),'polR']).ImageClickedFcn = {@ea_inc_polarity,handles,options,['k',num2str(k),'R']};
end

for k=1:S.numContacts
    if S.(['Ls',num2str(Lactive)]).(['k',num2str(k)]).pol==0 % off
        icon=fullfile(options.earoot, 'icons', ['empty',num2str(Lactive),'.png']);
    elseif S.(['Ls',num2str(Lactive)]).(['k',num2str(k)]).pol==1 % negative S1
        icon=fullfile(options.earoot, 'icons', ['minus',num2str(Lactive),'.png']);
    elseif S.(['Ls',num2str(Lactive)]).(['k',num2str(k)]).pol==2 % positive S1
        icon=fullfile(options.earoot, 'icons', ['plus',num2str(Lactive),'.png']);
    end
    handles.(['k',num2str(k),'polL']).ImageSource = icon;
    handles.(['k',num2str(k),'polL']).ImageClickedFcn = {@ea_inc_polarity,handles,options,['k',num2str(k),'L']};
end

for k=S.numContacts+1:16
    handles.(['k',num2str(k), 'polR']).Visible = 'off';
    handles.(['k',num2str(k), 'polL']).Visible = 'off';
end

% right case:
if S.(['Rs',num2str(Ractive)]).case.pol==0 % off
    icon=fullfile(options.earoot, 'icons', ['empty',num2str(Ractive),'.png']);
elseif S.(['Rs',num2str(Ractive)]).case.pol==1 % negative
    icon=fullfile(options.earoot, 'icons', ['minus',num2str(Ractive),'.png']);
elseif S.(['Rs',num2str(Ractive)]).case.pol==2 % positive
    icon=fullfile(options.earoot, 'icons', ['plus',num2str(Ractive),'.png']);
end
handles.CpolR.ImageSource = icon;
handles.CpolR.ImageClickedFcn = {@ea_inc_polarity,handles,options,'CpolR'};

% left case:
if S.(['Ls',num2str(Lactive)]).case.pol==0 % off
    icon=fullfile(options.earoot, 'icons', ['empty',num2str(Lactive),'.png']);
elseif S.(['Ls',num2str(Lactive)]).case.pol==1 % negative
    icon=fullfile(options.earoot, 'icons', ['minus',num2str(Lactive),'.png']);
elseif S.(['Ls',num2str(Lactive)]).case.pol==2 % positive
    icon=fullfile(options.earoot, 'icons', ['plus',num2str(Lactive),'.png']);
end
handles.CpolL.ImageSource = icon;
handles.CpolL.ImageClickedFcn = {@ea_inc_polarity,handles,options,'CpolL'};

%% add label

% set(handles.stimlabel,'String',S.label);

% check consistency with chosen VAT model.
% check consistency with chosen electrode model.
if ~isfield(options,'elspec')
    if isempty(elstruct(actpt).elmodel)
        error('Model is empty. Was the electrode segmentation fully executed in single patient mode?')
    end
    toptions=ea_resolve_elspec(elstruct(actpt));
    try
        options.elspec=toptions.elspec;
    catch
        keyboard
    end
end

ea_toggle_contacts(handles, S.numContacts);

%if strcmp(options.elspec.matfname,'boston_vercise_directed')
%    ea_error('VTA modeling for directed leads is not yet supported.');
%end

if get(handles.(['Rs',num2str(Ractive),'va']),'Value')==1 % Volt
    ea_show_percent(handles,1,'off'); % right hemisphere
else % Ampere
    ea_show_percent(handles,1,'on'); % right hemisphere
end
if get(handles.(['Ls',num2str(Ractive),'va']),'Value')==1 % Volt
    ea_show_percent(handles,2,'off'); % left hemisphere
else % Ampere
    ea_show_percent(handles,2,'on'); % left hemisphere
end

%% enable/disable panel based on sides that are present
is_side_present=arrayfun(@(xside) ~ea_arenopoints4side(elstruct(actpt).trajectory, xside), [1,2]);%First element is R, second is L
tabSide1 = findobj(handles.stimfig.Children, 'Type', 'uitab', 'Title', 'Right Hemisphere');
if is_side_present(1)>0 % check if R side is present
    set(findall(tabSide1, '-property', 'enable'), 'enable', 'on')
    %fix color (ensure they are reloaded correctly)
    handles.Rs1am.BackgroundColor=[1 1 1];%force redraw color
    handles.Rs1am.BackgroundColor=[0.953 0.871 0.733];%orange
    handles.Rs2am.BackgroundColor=[1 1 1];%force redraw color
    handles.Rs2am.BackgroundColor=[0.729 0.831 0.957];%blue
    handles.Rs3am.BackgroundColor=[1 1 1];%force redraw color
    handles.Rs3am.BackgroundColor=[0.925 0.839 0.839];%red
    handles.Rs4am.BackgroundColor=[1 1 1];%force redraw color
    handles.Rs4am.BackgroundColor=[0.757 0.867 0.776];%green
else
    set(findall(tabSide1, '-property', 'enable'), 'enable', 'off')
end

tabSide2 = findobj(handles.stimfig.Children, 'Type', 'uitab', 'Title', 'Left Hemisphere');
if is_side_present(2)>0 % check if L side is present
    set(findall(tabSide2, '-property', 'enable'), 'enable', 'on')
    %fix color (ensure they are reloaded correctly)
    handles.Ls1am.BackgroundColor=[1 1 1];%force redraw color
    handles.Ls1am.BackgroundColor=[0.953 0.871 0.733];%orange
    handles.Ls2am.BackgroundColor=[1 1 1];%force redraw color
    handles.Ls2am.BackgroundColor=[0.729 0.831 0.957];%blue
    handles.Ls3am.BackgroundColor=[1 1 1];%force redraw color
    handles.Ls3am.BackgroundColor=[0.925 0.839 0.839];%red
    handles.Ls4am.BackgroundColor=[1 1 1];%force redraw color
    handles.Ls4am.BackgroundColor=[0.757 0.867 0.776];%green
else
    set(findall(tabSide2, '-property', 'enable'), 'enable', 'off')
end

% if contains(model, 'OSS-DBS') && vatsettings.butenko_calcPAM
%     handles.Rs2am.Enable = "off";
%     handles.Rs3am.Enable = "off";
%     handles.Rs4am.Enable = "off";
%     handles.Ls2am.Enable = "off";
%     handles.Ls3am.Enable = "off";
%     handles.Ls4am.Enable = "off";
% else
%     handles.Rs2am.Enable = "on";
%     handles.Rs3am.Enable = "on";
%     handles.Rs4am.Enable = "on";
%     handles.Ls2am.Enable = "on";
%     handles.Ls3am.Enable = "on";
%     handles.Ls4am.Enable = "on";
% end

switch model
    case 'SimBio/FieldTrip (see Horn 2017)'
        ea_hide_impedance(handles);
        set(handles.estimateInTemplate,'Visible','on');
        S.monopolarmodel=0;
        ea_enable_vas(handles,options);
        ea_toggle_pulsewidth(handles, 'off');
        set(handles.betawarning,'visible','on');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
    case 'Maedler 2012'
        ea_show_impedance(handles);
        ea_setprefs('vatsettings.estimateInTemplate', 1);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=1;
        ea_disable_vas(handles,options);
        ea_toggle_pulsewidth(handles, 'off');
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','off');
        set(handles.addStimSet,'visible','off');
    case 'Kuncel 2008'
        ea_hide_impedance(handles);
        ea_setprefs('vatsettings.estimateInTemplate', 1);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=1;
        ea_disable_vas(handles,options);
        ea_toggle_pulsewidth(handles, 'off');
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','off');
        set(handles.addStimSet,'visible','off');
    case 'Dembek 2017'
        ea_show_impedance(handles);
        ea_setprefs('vatsettings.estimateInTemplate', 1);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=1;
        ea_enable_vas(handles,options);
        ea_toggle_pulsewidth(handles, 'off');
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
    case 'Fastfield (Baniasadi 2020)'
        ea_show_impedance(handles);
        ea_setprefs('vatsettings.estimateInTemplate', 1);
        set(handles.estimateInTemplate,'Visible','off');
        S.monopolarmodel=0;
        ea_enable_vas(handles,options);
        ea_toggle_pulsewidth(handles, 'off');
        set(handles.betawarning,'visible','off');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
    case 'OSS-DBS (Butenko 2020)'
        ea_hide_impedance(handles);
        set(handles.estimateInTemplate,'Visible','on');
        S.monopolarmodel=0;

        ea_enable_vas(handles,options);
        ea_toggle_pulsewidth(handles, 'on');

        set(handles.betawarning,'visible','on');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');
end

% Set estimate in template to 0 when visualize in native space
if options.native
    ea_setprefs('vatsettings.estimateInTemplate', 0);
    set(handles.estimateInTemplate, 'Visible', 'off');
end

S.model=model;

ea_savestimulation(S,options);
setappdata(handles.stimfig,'S',S);


function ea_show_percent(handles,side,onoff)
S = getappdata(handles.stimfig, 'S');
switch side
    case 1
        sidestr='R';
    case 2
        sidestr='L';
end

for k=1:S.numContacts
    set(handles.(['k',num2str(k),'u',sidestr]), 'visible', onoff);
end

set(handles.(['Cu', sidestr]),'visible',onoff);

for i = 1:ceil(S.numContacts/4)
    set(handles.(['perctext', num2str(i), sidestr]), 'visible', onoff);
end


function ea_toggle_contacts(handles, numel)
for k=1:numel
    set(handles.(['k',num2str(k),'uR']), 'visible', 'on');
    set(handles.(['k',num2str(k),'uL']), 'visible', 'on');
    set(handles.(['k',num2str(k),'imR']), 'visible', 'on');
    set(handles.(['k',num2str(k),'imL']), 'visible', 'on');
    set(handles.(['k',num2str(k),'txtR']),'visible', 'on');
    set(handles.(['k',num2str(k),'txtL']),'visible', 'on');
    set(handles.(['k',num2str(k),'polR']),'visible', 'on');
    set(handles.(['k',num2str(k),'polL']),'visible', 'on');
end

for k=numel+1:16
    set(handles.(['k',num2str(k),'uR']), 'visible', 'off');
    set(handles.(['k',num2str(k),'uL']), 'visible', 'off');
    set(handles.(['k',num2str(k),'imR']), 'visible', 'off');
    set(handles.(['k',num2str(k),'imL']), 'visible', 'off');
    set(handles.(['k',num2str(k),'txtR']),'visible', 'off');
    set(handles.(['k',num2str(k),'txtL']),'visible', 'off');
    set(handles.(['k',num2str(k),'polR']),'visible', 'off');
    set(handles.(['k',num2str(k),'polL']),'visible', 'off');
end

for i = 1:ceil(numel/4)
    set(handles.(['perctext', num2str(i), 'R']), 'visible', 'on');
    set(handles.(['perctext', num2str(i), 'L']), 'visible', 'on');
    set(handles.(['kohmtext', num2str(i), 'R']),'visible', 'on');
    set(handles.(['kohmtext', num2str(i), 'L']),'visible', 'on');
end

for i = ceil(numel/4)+1:4
    set(handles.(['perctext', num2str(i), 'R']), 'visible', 'off');
    set(handles.(['perctext', num2str(i), 'L']), 'visible', 'off');
    set(handles.(['kohmtext', num2str(i), 'R']),'visible', 'off');
    set(handles.(['kohmtext', num2str(i), 'L']),'visible', 'off');
end


function ea_disable_vas(handles,options)
RL={'R','L'};
for iside=1:length(options.sides)
    side=options.sides(iside);

    for Rva=1:4
        set(handles.([RL{side},'s',num2str(Rva),'va']),'enable','off');
        set(handles.([RL{side},'s',num2str(Rva),'va']),'value',1);
    end
end


function ea_enable_vas(handles,options)
RL={'R','L'};
for iside=1:length(options.sides)
    side=options.sides(iside);

    for Rva=1:4
        set(handles.([RL{side},'s',num2str(Rva),'va']),'enable','on');
    end
end


function ea_toggle_pulsewidth(handles, status)
handles.pulseWidthLabel_R.Visible = status;
handles.pulseWidthTextbox_R.Visible = status;
handles.usLabel_R.Visible = status;
handles.pulseWidthLabel_L.Visible = status;
handles.pulseWidthTextbox_L.Visible = status;
handles.usLabel_L.Visible = status;


function ea_hide_impedance(handles)
for k=1:16
    set(handles.(['k', num2str(k), 'imR']), 'Visible', 'off');
    set(handles.(['k', num2str(k), 'imL']), 'Visible', 'off');
end

for ohm=1:4
    set(handles.(['kohmtext', num2str(ohm), 'R']), 'Visible', 'off');
    set(handles.(['kohmtext', num2str(ohm), 'L']), 'Visible', 'off');
end


function ea_show_impedance(handles)
S = getappdata(handles.stimfig, 'S');
for k=1:S.numContacts
    set(handles.(['k', num2str(k), 'imR']), 'Visible', 'on');
    set(handles.(['k', num2str(k), 'imL']), 'Visible', 'on');
end

for ohm=1:ceil(S.numContacts/4)
    set(handles.(['kohmtext', num2str(ohm), 'R']), 'Visible', 'on');
    set(handles.(['kohmtext', num2str(ohm), 'L']), 'Visible', 'on');
end


function S=ea_redistribute_voltage(S,changedobj)
conts = strcat('k', arrayfun(@num2str, 1:S.numContacts, 'Uni', 0));
contsCase = [conts, 'case'];
pulseWidthTextbox = {'pulseWidthTextbox_R', 'pulseWidthTextbox_L'};

if ischar(changedobj) % different polarity on the block
    if ismember(changedobj, pulseWidthTextbox)
        return;
    end

    if endsWith(changedobj, 'R')
        side = 1;
    elseif endsWith(changedobj, 'L')
        side = 2;
    end
    sidestr = changedobj(end);

    if startsWith(changedobj, 'C', 'IgnoreCase', true)
        changedobj = 'case';
    elseif startsWith(changedobj, 'k')
        changedobj = regexp(changedobj, '^k\d+', 'match', 'once');
    end

    % check polarity of changed object:
    polchanged=S.([sidestr,'s',num2str(S.active(side))]).(changedobj).pol;

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).pol=0;
            S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).perc=0;
        end
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).pol=1;
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).perc=100;
        return
    else
        % if S.([sidec,'s',num2str(S.active(side))]).va==2 % ampere only allows one anode and one cathode
        %     for c=1:length(contsCase)
        %         if S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol==polchanged % same polarity as changed object
        %             S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=ea_swappol(polchanged);
        %             S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=100;
        %         else
        %             S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=0;
        %             S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=0;
        %         end
        %     end
        %     S.([sidec,'s',num2str(S.active(side))]).(changedobj).pol=1;
        %     S.([sidec,'s',num2str(S.active(side))]).(changedobj).perc=100;
        % end
    end

    if polchanged==0
        % set changed contacts percentage to zero:
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).perc=0;
    else
        % determine how many other nodes with this polarity exist:
        divby=1;
        contacts={};
        for c=1:length(conts)
            if S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).pol==polchanged
                if ~strcmp(conts{c},changedobj)
                    %voltages{divby}=S.([sidec,'s',num2str(S.active(side))]).(conts{c}).perc;
                    contacts{divby}=conts{c};
                    divby=divby+1;
                end
            end
        end

        if S.([sidestr,'s',num2str(S.active(side))]).case.pol==polchanged
            if ~strcmp(changedobj,'case')
                contacts{divby}='case';
                divby=divby+1;
            end
        end
        % add case to caCuLlation.

        % set changed contacts percentage:
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).perc=100/divby;

        % reduce all other contacts percentages:

        try divby=divby/length(contacts); end
        for c=1:length(contacts)
            S.([sidestr,'s',num2str(S.active(side))]).(contacts{c}).perc=S.([sidestr,'s',num2str(S.active(side))]).(contacts{c}).perc/divby;
        end
    end

    % now clean up mess from polarity that the contact used to have..

    polchanged=ea_polminus(polchanged);
    sumpercs=0;

    if polchanged % polarization has changed from negative to positive. clean up negatives. or changed from positive to off. clean up positives.
        contacts={};
        cnt=0;
        for c=1:length(conts)
            if S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).pol==polchanged
                if ~strcmp(conts{c},changedobj)
                    %voltages{divby}=S.([sidec,'s',num2str(S.active(side))]).(Rconts{con}).perc;

                    cnt=cnt+1;
                    contacts{cnt}=conts{c};
                    sumpercs=sumpercs+S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).perc;
                end
            end
        end
        % add case to caCuLlation:
        if S.([sidestr,'s',num2str(S.active(side))]).case.pol==polchanged
            if ~strcmp(changedobj,'case')
                cnt=cnt+1;
                contacts{cnt}='case';
                sumpercs=sumpercs+S.([sidestr,'s',num2str(S.active(side))]).case.perc;
            end
        end

        multby=(100/sumpercs);
        if cnt
            for c=1:length(contacts)
                S.([sidestr,'s',num2str(S.active(side))]).(contacts{c}).perc=S.([sidestr,'s',num2str(S.active(side))]).(contacts{c}).perc*multby;
            end
        end
    end

else % voltage percentage changed
    changedobj=get(changedobj, 'Tag');
    if ismember(changedobj, pulseWidthTextbox)
        return;
    end

    if endsWith(changedobj, 'R')
        side = 1;
    elseif endsWith(changedobj, 'L')
        side = 2;
    end
    sidestr = changedobj(end);

    if startsWith(changedobj, 'C', 'IgnoreCase', true)
        changedobj = 'case';
    elseif startsWith(changedobj, 'k')
        changedobj = regexp(changedobj, '^k\d+', 'match', 'once');
    end

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).pol=0;
            S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).perc=0;
        end
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).pol=1;
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).perc=100;

        return
    end

    % check polarity of changed object:
    try
        polchanged=S.([sidestr,'s',num2str(S.active(side))]).(changedobj).pol;
    catch
        keyboard
    end

    if polchanged==0 % set changed contacts polarity to negative
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).pol=1;
        polchanged=1;
    end

    % determine how many other nodes with this polarity exist:
    divby=1;
    contacts={};
    sumpercent=0;
    for c=1:length(conts)
        if S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).pol==polchanged
            if ~strcmp(conts{c},changedobj)
                sumpercent=sumpercent+S.([sidestr,'s',num2str(S.active(side))]).(conts{c}).perc;
                contacts{divby}=conts{c};
                divby=divby+1;
            end
        end
    end

    % add case to caCuLlation.
    if S.([sidestr,'s',num2str(S.active(side))]).case.pol==polchanged
        if ~strcmp(changedobj,'case')
            contacts{divby}='case';
            divby=divby+1;
        end
    end

    if divby==1 % only one contact -> set to 100 percent.
        S.([sidestr,'s',num2str(S.active(side))]).(changedobj).perc=100;
    end

    % reduce all other contacts percentages:
    divby=sumpercent/(100-S.([sidestr,'s',num2str(S.active(side))]).(changedobj).perc);

    for c=1:length(contacts)
        S.([sidestr,'s',num2str(S.active(side))]).(contacts{c}).perc=S.([sidestr,'s',num2str(S.active(side))]).(contacts{c}).perc/divby;
    end
end


function opol=ea_polminus(pol)
if pol==1
    opol=0;
elseif pol==2
    opol=1;
elseif pol==0
    opol=2;
end


function opol=ea_swappol(pol)
if pol==1
    opol=2;
elseif pol==2
    opol=1;
else
    opol=0;
end


function ea_inc_polarity(~,~,handles,options,ID)

S=getappdata(handles.stimfig,'S');

if endsWith(ID, 'R')
    side = 1;
elseif endsWith(ID, 'L')
    side = 2;
end
sidestr = ID(end);

if startsWith(ID, 'C', 'IgnoreCase', true)
    gID = 'case';
elseif startsWith(ID, 'k')
    gID = ID(1:end-1);
end

cycles=[0,1,2];

try
    oldval=S.([sidestr,'s',num2str(S.active(side))]).(gID).pol;
catch
    keyboard
end
[~,newval]=ismember(oldval,cycles);
newval=newval+1;
if newval>length(cycles)
    newval=1;
end

newval=cycles(newval);

S.([sidestr,'s',num2str(S.active(side))]).(gID).pol=newval;

% now check if any other contact is left with the old polarity

anycontactpositive=0;
anycontactnegative=0;
for k=1:S.numContacts
    if S.([sidestr,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==1
        anycontactnegative=1;
    elseif S.([sidestr,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==2
        anycontactpositive=1;
    end
end

% also check case
if S.([sidestr,'s',num2str(S.active(side))]).case.pol==1
    anycontactnegative=1;
elseif S.([sidestr,'s',num2str(S.active(side))]).case.pol==2
    anycontactpositive=1;
end

if anycontactnegative && anycontactpositive % only then save results..
    setappdata(handles.stimfig,'S',S);
end
ea_refreshguisp(handles,options,ID);
