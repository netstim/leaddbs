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
    S=ea_initializeS(stimlabel,options,handles);
    try S.model=gSv.vatmodel; end
    setappdata(handles.stimfig,'stimlabel',S.label);
else
    if isempty(S.Rs1)
        S=ea_initializeS(stimlabel,options,handles);
        try S.model=gSv.vatmodel; end
        setappdata(handles.stimfig,'stimlabel',S.label);
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
if contains(model, 'OSS-DBS') && vatsettings.butenko_calcPAM
    handles.Rs2am.Value = 0;
    handles.Rs3am.Value = 0;
    handles.Rs4am.Value = 0;
    handles.Ls2am.Value = 0;
    handles.Ls3am.Value = 0;
    handles.Ls4am.Value = 0;
    S.active = [1 1];
end

Ractive=S.active(1);
Lactive=S.active(2);

if nargin==3
    if ischar(varargin{3})
        switch varargin{3}
            case {'Rcase'}
                ks={'k0','k1','k2','k3','k4','k5','k6','k7'};
                sidec='R'; side=1;
                S=ea_redistribute_voltage(S,varargin{3});
                if S.(['Rs',num2str(Ractive)]).case.pol==1

                    S.(['Rs',num2str(Ractive)]).case.pol=2;
                    S=ea_redistribute_voltage(S,varargin{3});
                    for k=0:7
                        if S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==2
                            S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol=0;
                            S=ea_redistribute_voltage(S,['k',num2str(k)]);
                        end
                    end
                end
            case {'Lcase'}
                ks={'k8','k9','k10','k11','k12','k13','k14','k15'};
                sidec='L'; side=2;
                S=ea_redistribute_voltage(S,varargin{3});
                if S.(['Ls',num2str(Lactive)]).case.pol==1
                    S.(['Ls',num2str(Lactive)]).case.pol=2;
                    S=ea_redistribute_voltage(S,varargin{3});
                    for k=8:15
                        if S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol==2
                            S.([sidec,'s',num2str(S.active(side))]).(['k',num2str(k)]).pol=0;
                            S=ea_redistribute_voltage(S,['k',num2str(k)]);
                        end
                    end
                end
            otherwise
                S=ea_redistribute_voltage(S,varargin{3});

                switch varargin{3}
                    case {'k0','k1','k2','k3','k4','k5','k6','k7'}
                        sidec='R'; side=1;
                    case {'k8','k9','k10','k11','k12','k13','k14','k15'}
                        sidec='L'; side=2;
                end
                if S.([sidec,'s',num2str(S.active(side))]).(varargin{3}).pol==2 && S.([sidec,'s',num2str(S.active(side))]).case.pol==2
                    S.([sidec,'s',num2str(S.active(side))]).case.pol=0;
                    S=ea_redistribute_voltage(S,[sidec,'case']);
                end
        end
    else
        S=ea_redistribute_voltage(S,varargin{3});
    end
end

if contains(model, 'OSS-DBS')
    vatsettings.butenko_pulseWidth = ea_getprefs('vatsettings.butenko_pulseWidth');
    for source = 1:4
        if eval(['~isfield(', 'S.Rs', num2str(source), ',''pulseWidth'')'])
            eval(['S.Rs',num2str(source),'.pulseWidth=',num2str(vatsettings.butenko_pulseWidth),';']);
        end
        if eval(['~isfield(', 'S.Ls', num2str(source), ',''pulseWidth'')'])
            eval(['S.Ls',num2str(source),'.pulseWidth=',num2str(vatsettings.butenko_pulseWidth),';']);
        end
    end
else
    for source = 1:4
        if eval(['isfield(', 'S.Rs', num2str(source), ',''pulseWidth'')'])
            eval(['S.Rs', num2str(source), ' = rmfield(', 'S.Rs', num2str(source), ',''pulseWidth'');']);
        end
        if eval(['isfield(', 'S.Ls', num2str(source), ',''pulseWidth'')'])
            eval(['S.Ls', num2str(source), ' = rmfield(', 'S.Ls', num2str(source), ',''pulseWidth'');']);
        end
    end
end

setappdata(handles.stimfig,'S',S);

% set stim amplitudes
for source=1:4
    S.amplitude{1}(source)=S.(['Rs',num2str(source)]).amp;

    set(eval(['handles.Rs',num2str(source),'am']),'String',num2str(S.amplitude{1}(source)));
    set(eval(['handles.Rs',num2str(source),'va']),'Value',eval(['S.Rs',num2str(source),'.va']));

    %if eval(['S.Rs',num2str(source),'.amp']) % check if a valid +/- combination is active, if not set defaults.
    anycontactpositive=0; anycontactnegative=0;
    for k=0:7
        if eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol==1'])
            anycontactnegative=1;
        elseif eval(['S.Rs',num2str(source),'.k',num2str(k),'.pol==2'])
            anycontactpositive=2;
        end
    end

    if ~anycontactnegative
        eval(['S.Rs',num2str(source),'.k1.pol=1;']);
        eval(['S.Rs',num2str(source),'.k1.perc=100;']);
    end

    if ~anycontactpositive
        eval(['S.Rs',num2str(source),'.case.pol=2;']);
        eval(['S.Rs',num2str(source),'.case.perc=100;']);
    end
    %end
end

for source=1:4
    S.amplitude{2}(source)=S.(['Ls',num2str(source)]).amp;

    set(eval(['handles.Ls',num2str(source),'am']),'String',num2str(S.amplitude{2}(source)));
    set(eval(['handles.Ls',num2str(source),'va']),'Value',eval(['S.Ls',num2str(source),'.va']));

    % if eval(['S.Ls',num2str(source),'.amp']) % check if a valid +/- combination is active, if not set defaults.
    anycontactpositive=0; anycontactnegative=0;
    for k=8:15
        if eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol==1'])
            anycontactnegative=1;
        elseif eval(['S.Ls',num2str(source),'.k',num2str(k),'.pol==2'])
            anycontactpositive=1;
        end
    end

    if ~anycontactnegative
        eval(['S.Ls',num2str(source),'.k9.pol=1;']);
        eval(['S.Ls',num2str(source),'.k9.perc=100;']);
    end
    if ~anycontactpositive
        eval(['S.Ls',num2str(source),'.case.pol=2;']);
        eval(['S.Ls',num2str(source),'.case.perc=100;']);
    end
    % end
end

%% model to handles: all GUI elements.
source=Ractive;
for k=0:7
    val=eval(['S.Rs',num2str(source),'.k',num2str(k),'.perc']);
    set(eval(['handles.k',num2str(k),'u']),'String',num2str(val));

    val=eval(['S.Rs',num2str(source),'.k',num2str(k),'.imp']);
    set(eval(['handles.k',num2str(k),'im']),'String',num2str(val));
end

% set case
set(handles.RCu,'String',num2str(eval(['S.Rs',num2str(source),'.case.perc'])));

source=Lactive;
for k=8:15
    val=eval(['S.Ls',num2str(source),'.k',num2str(k),'.perc']);
    set(eval(['handles.k',num2str(k),'u']),'String',num2str(val));

    val=eval(['S.Ls',num2str(source),'.k',num2str(k),'.imp']);
    set(eval(['handles.k',num2str(k),'im']),'String',num2str(val));
end

% set case
set(handles.LCu,'String',num2str(eval(['S.Ls',num2str(source),'.case.perc'])));

if contains(model, 'OSS-DBS')
    handles.pulseWidthTextbox_R.String = num2str(eval(['S.Rs',num2str(S.active(1)),'.pulseWidth']));
    handles.pulseWidthTextbox_L.String = num2str(eval(['S.Ls',num2str(S.active(2)),'.pulseWidth']));
end

%% model to handles: Axes objects:
for k=0:7
    if eval(['S.Rs',num2str(Ractive),'.k',num2str(k),'.pol==0']) % off
        im=ea_get_icn(['empty',num2str(Ractive)]);
    elseif eval(['S.Rs',num2str(Ractive),'.k',num2str(k),'.pol==1']) % negative S1
        im=ea_get_icn(['minus',num2str(Ractive)]);
    elseif eval(['S.Rs',num2str(Ractive),'.k',num2str(k),'.pol==2']) % positive S1
        im=ea_get_icn(['plus',num2str(Ractive)]);
    end
    set(0,'CurrentFigure',handles.stimfig);
    set(handles.stimfig,'CurrentAxes',eval(['handles.k',num2str(k),'ax']));
    h=image(im);
    set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,['k',num2str(k)]});
    axis off;
    axis equal;
end

for k=8:15
    if eval(['S.Ls',num2str(Lactive),'.k',num2str(k),'.pol==0']) % off
        im=ea_get_icn(['empty',num2str(Lactive)]);
    elseif eval(['S.Ls',num2str(Lactive),'.k',num2str(k),'.pol==1']) % negative S1
        im=ea_get_icn(['minus',num2str(Lactive)]);
    elseif eval(['S.Ls',num2str(Lactive),'.k',num2str(k),'.pol==2']) % positive S1
        im=ea_get_icn(['plus',num2str(Lactive)]);
    end
    set(0,'CurrentFigure',handles.stimfig);
    set(handles.stimfig,'CurrentAxes',eval(['handles.k',num2str(k),'ax']));
    h=image(im);
    set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,['k',num2str(k)]});
    axis off;
    axis equal;
end

% right case:
if eval(['S.Rs',num2str(Ractive),'.case.pol==0']) % off
    im=ea_get_icn(['empty',num2str(Ractive)]);
elseif eval(['S.Rs',num2str(Ractive),'.case.pol==1']) % negative
    im=ea_get_icn(['minus',num2str(Ractive)]);
elseif eval(['S.Rs',num2str(Ractive),'.case.pol==2']) % positive
    im=ea_get_icn(['plus',num2str(Ractive)]);
end
set(0,'CurrentFigure',handles.stimfig);
set(handles.stimfig,'CurrentAxes',handles.RCax);

h=image(im);
set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,'Rcase'});
axis off;
axis equal;

% left case:
if eval(['S.Ls',num2str(Lactive),'.case.pol==0']) % off
    im=ea_get_icn(['empty',num2str(Lactive)]);
elseif eval(['S.Ls',num2str(Lactive),'.case.pol==1']) % negative
    im=ea_get_icn(['minus',num2str(Lactive)]);
elseif eval(['S.Ls',num2str(Lactive),'.case.pol==2']) % positive
    im=ea_get_icn(['plus',num2str(Lactive)]);
end
set(0,'CurrentFigure',handles.stimfig);
set(handles.stimfig,'CurrentAxes',handles.LCax);

h=image(im);
set(h,'ButtonDownFcn',{@ea_inc_polarity,handles,options,'Lcase'});
axis off;
axis equal;

%% add label

%set(handles.stimlabel,'String',S.label);

%% check consistency with chosen VAT model.
%% check consistency with chosen electrode model.
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

if options.elspec.numel > 8
    warning('Only electrode with less than 8 contacts are fully supported.');
else
    ea_toggle_contacts(handles, options.elspec.numel);
end

%if strcmp(options.elspec.matfname,'boston_vercise_directed')
%    ea_error('VTA modeling for directed leads is not yet supported.');
%end

if get(handles.(['Rs',num2str(Ractive),'va']),'Value')==1 % Volt
    ea_show_percent(handles,options,1,'off'); % right hemisphere
else % Ampere
    ea_show_percent(handles,options,1,'on'); % right hemisphere
end
if get(handles.(['Ls',num2str(Ractive),'va']),'Value')==1 % Volt
    ea_show_percent(handles,options,2,'off'); % left hemisphere
else % Ampere
    ea_show_percent(handles,options,2,'on'); % left hemisphere
end

%% enable/disable panel based on sides that are present
is_side_present=arrayfun(@(xside) ~ea_arenopoints4side(elstruct(actpt).trajectory, xside), [1,2]);%First element is R, second is L
if is_side_present(1)>0%check if R side is present
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'on')
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
    set(findall(handles.uipanel2, '-property', 'enable'), 'enable', 'off')
end
if is_side_present(2)>0%check if L side is present
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'on')
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
    set(findall(handles.uipanel3, '-property', 'enable'), 'enable', 'off')
end

if contains(model, 'OSS-DBS') && vatsettings.butenko_calcPAM
    handles.Rs2am.Enable = "off";
    handles.Rs3am.Enable = "off";
    handles.Rs4am.Enable = "off";
    handles.Ls2am.Enable = "off";
    handles.Ls3am.Enable = "off";
    handles.Ls4am.Enable = "off";
else
    handles.Rs2am.Enable = "on";
    handles.Rs3am.Enable = "on";
    handles.Rs4am.Enable = "on";
    handles.Ls2am.Enable = "on";
    handles.Ls3am.Enable = "on";
    handles.Ls4am.Enable = "on";
end

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
        if vatsettings.butenko_calcPAM
            Rs1va = handles.Rs1va.Value;
            Ls1va = handles.Ls1va.Value;
            ea_disable_vas(handles,options);
            handles.Rs1va.Enable = "on";
            handles.Rs1va.Value = Rs1va;
            handles.Ls1va.Enable = "on";
            handles.Ls1va.Value = Ls1va;
            ea_toggle_pulsewidth(handles, 'on');
        else
            ea_enable_vas(handles,options);
            ea_toggle_pulsewidth(handles, 'on');
        end
        set(handles.betawarning,'visible','on');
        set(handles.settings,'visible','on');
        set(handles.addStimSet,'visible','off');

end
S.model=model;

ea_savestimulation(S,options);
setappdata(handles.stimfig,'S',S);


function ea_show_percent(handles,options,side,onoff)

switch side
    case 1
        sel=0:7;
        sidestr='R';
        ptval=1;
    case 2
        sel=8:15;
        sidestr='L';
        ptval=3;
end

% Only support up to 8 contacts for now
if options.elspec.numel<=8
    sel=sel(1:options.elspec.numel);
else
    sel=sel(1:8);
end

for k=sel
    set(handles.(['k',num2str(k),'u']),'visible',onoff);
end

set(handles.([sidestr,'Cu']),'visible',onoff);

set(handles.(['perctext',num2str(ptval)]),'visible',onoff);
if options.elspec.numel>4
    set(handles.(['perctext',num2str(ptval+1)]),'visible',onoff);
end


function ea_toggle_contacts(handles, numel)

if numel == 8
    status = 'on';
else
    status = 'off';
end


for k=[numel:7,numel+8:15]
    set(handles.(['k',num2str(k),'u']),'visible',status);
    set(handles.(['k',num2str(k),'im']),'visible',status);
    set(handles.(['k',num2str(k),'txt']),'visible',status);
    handles2hide = [get(handles.(['k',num2str(k),'ax']),'Children')];
    set(handles2hide,'visible',status)
end

set(handles.perctext2,'visible',status);
set(handles.kohmtext2,'visible',status);
set(handles.perctext4,'visible',status);
set(handles.kohmtext4,'visible',status);


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

for k=0:15
    eval(['set(handles.k',num2str(k),'im,''visible'',''off'');']);
end

for ohm=1:4
    eval(['set(handles.kohmtext',num2str(ohm),',''visible'',''off'');']);
end


function ea_show_impedance(handles)

for k=0:15
    eval(['set(handles.k',num2str(k),'im,''visible'',''on'');']);
end

for ohm=1:4
    eval(['set(handles.kohmtext',num2str(ohm),',''visible'',''on'');']);
end


function S=ea_redistribute_voltage(S,changedobj)
Rconts={'k0','k1','k2','k3','k4','k5','k6','k7'};
Lconts={'k8','k9','k10','k11','k12','k13','k14','k15'};
LcontsCase=[Lconts,{'case'}];
RcontsCase=[Rconts,{'case'}];
pulseWidthTextbox = {'pulseWidthTextbox_R', 'pulseWidthTextbox_L'};
if ischar(changedobj) % different polarity on the block
    switch changedobj
        case Rconts
            conts=Rconts;
            contsCase=RcontsCase;
            sidec='R';
            side=1;
        case Lconts
            conts=Lconts;
            contsCase=LcontsCase;
            sidec='L';
            side=2;
        case 'Rcase'
            conts=Rconts;
            changedobj='case';
            contsCase=RcontsCase;

            side=1;
            sidec='R';
        case 'Lcase'
            conts=Lconts;
            contsCase=LcontsCase;

            changedobj='case';
            side=2;
            sidec='L';
        case pulseWidthTextbox
            return;
    end

    % check polarity of changed object:
    polchanged=eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol']);

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.pol=0;']);
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.perc=0;']);
        end
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);

        return

    else
        %         if S.([sidec,'s',num2str(S.active(side))]).va==2 % ampere only allows one anode and one cathode
        %             for c=1:length(contsCase)
        %
        %                 if S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol==polchanged % same polarity as changed object
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=ea_swappol(polchanged);
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=100;
        %                 else
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).pol=0;
        %                     S.([sidec,'s',num2str(S.active(side))]).(contsCase{c}).perc=0;
        %                 end
        %             end
        %             S.([sidec,'s',num2str(S.active(side))]).(changedobj).pol=1;
        %             S.([sidec,'s',num2str(S.active(side))]).(changedobj).perc=100;
        %         end
    end

    if polchanged==0
        % set changed contacts percentage to zero:
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=0;']);
    else
        % determine how many other nodes with this polarity exist:
        divby=1;
        contacts={};
        for con=1:length(conts)
            if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
                if ~strcmp(conts{con},changedobj)
                    %voltages{divby}=eval(['S.Rs',num2str(S.active(side)),'.',Rconts{con},'.perc']);
                    contacts{divby}=conts{con};
                    divby=divby+1;
                end
            end
        end

        if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
            if ~strcmp(changedobj,'case')
                contacts{divby}='case';
                divby=divby+1;
            end
        end
        % add case to calculation.

        % set changed contacts percentage:
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100/divby;']);

        % reduce all other contacts percentages:

        try divby=divby/length(contacts); end
        for c=1:length(contacts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
                'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc/divby;']);
        end
    end

    % now clean up mess from polarity that the contact used to have..

    polchanged=ea_polminus(polchanged);
    sumpercs=0;

    if polchanged % polarization has changed from negative to positive. clean up negatives. or changed from positive to off. clean up positives.
        contacts={};
        cnt=0;
        for con=1:length(conts)
            if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
                if ~strcmp(conts{con},changedobj)
                    %voltages{divby}=eval(['S.Rs',num2str(S.active(side)),'.',Rconts{con},'.perc']);

                    cnt=cnt+1;
                    contacts{cnt}=conts{con};
                    sumpercs=sumpercs+eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.perc']);
                end
            end
        end
        % add case to calculation:
        if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
            if ~strcmp(changedobj,'case')
                cnt=cnt+1;
                contacts{cnt}='case';
                sumpercs=sumpercs+eval(['S.',sidec,'s',num2str(S.active(side)),'.case.perc']);
            end
        end

        multby=(100/sumpercs);
        if cnt
            for c=1:length(contacts)
                eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
                    'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc*multby;']);
            end
        end
    end

else % voltage percentage changed
    changedobj=get(changedobj,'Tag');
    if ismember(changedobj, {'pulseWidthTextbox_R', 'pulseWidthTextbox_L'})
        return;
    end
    changedobj=changedobj(1:end-1);

    switch changedobj
        case Rconts
            conts=Rconts;
            sidec='R';
            side=1;
        case Lconts
            conts=Lconts;
            sidec='L';
            side=2;
        case 'RC'
            conts=Rconts;
            changedobj='case';
            side=1;
            sidec='R';
        case 'LC'
            conts=Lconts;
            changedobj='case';
            side=2;
            sidec='L';
    end

    % check for monopolar models:
    if S.monopolarmodel % these allow only 1 active anode contact per model.
        for c=1:length(conts)
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.pol=0;']);
            eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{c},'.perc=0;']);
        end
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);

        return
    end

    % check polarity of changed object:
    try
        polchanged=eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol']);
    catch
        keyboard
    end

    if polchanged==0 % set changed contacts polarity to negative
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.pol=1;']);
        polchanged=1;
    end

    % determine how many other nodes with this polarity exist:
    divby=1;
    contacts={};
    sumpercent=0;
    for con=1:length(conts)
        if eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.pol==polchanged'])
            if ~strcmp(conts{con},changedobj)
                sumpercent=sumpercent+eval(['S.',sidec,'s',num2str(S.active(side)),'.',conts{con},'.perc']);
                contacts{divby}=conts{con};
                divby=divby+1;
            end
        end
    end

    % add case to calculation.
    if eval(['S.',sidec,'s',num2str(S.active(side)),'.case.pol==polchanged'])
        if ~strcmp(changedobj,'case')
            contacts{divby}='case';
            divby=divby+1;
        end
    end

    if divby==1 % only one contact -> set to 100 percent.
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc=100;']);
    end

    % reduce all other contacts percentages:
    divby=sumpercent/(100-eval(['S.',sidec,'s',num2str(S.active(side)),'.',changedobj,'.perc']));

    for c=1:length(contacts)
        eval(['S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc=',...
            'S.',sidec,'s',num2str(S.active(side)),'.',contacts{c},'.perc/divby;']);
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


function ea_inc_polarity(h,h2,handles,options,ID)

S=getappdata(handles.stimfig,'S');

switch ID
    case {'k0','k1','k2','k3','k4','k5','k6','k7'}
        side=1;
        sidec='R';
        gID=ID;
        ks=0:7;
    case {'k8','k9','k10','k11','k12','k13','k14','k15'}
        side=2;
        sidec='L';
        gID=ID;
        ks=8:15;
    case 'Rcase'
        gID='case';
        side=1;
        sidec='R';
        ks=0:7;
    case 'Lcase'
        gID='case';
        side=2;
        sidec='L';
        ks=8:15;
end

cycles=[0,1,2];

try
    oldval=eval(['S.',sidec,'s',num2str(S.active(side)),'.',gID,'.pol']);
catch
    keyboard
end
[~,newval]=ismember(oldval,cycles);
newval=newval+1;
if newval>length(cycles)
    newval=1;
end

newval=cycles(newval);

eval(['S.',sidec,'s',num2str(S.active(side)),'.',gID,'.pol=',num2str(newval),';']);

% now check if any other contact is left with the old polarity

anycontactpositive=0; anycontactnegative=0;
for k=ks
    if eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.k',num2str(k),'.pol==1'])
        anycontactnegative=1;
    elseif eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.k',num2str(k),'.pol==2'])
        anycontactpositive=1;
    end
end

% also check case
if eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.case.pol==1'])
    anycontactnegative=1;
elseif eval(['S.',sidec,'s',num2str(num2str(S.active(side))),'.case.pol==2'])
    anycontactpositive=1;
end

if anycontactnegative && anycontactpositive % only then save results..
    setappdata(handles.stimfig,'S',S);
end
ea_refreshguisp(handles,options,ID);
