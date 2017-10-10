function ea_manualreconstruction(mcfig,patientname,options)
% This is the function that enables the user to manually correct the
% electrode positions that have been reconstructed by the algorithm. A
% handle for the figure, the coordinates of the contacts in millimeter
% notation, optionally manually measured coordinates, the fitted line in
% form of a 1x2 cell each containing a nx3 matrix that describes the line,
% the full path to the coronar nifti file, the name of the patient and the
% usual options struct must be handed to the function as parameters.
%
% Output parameters are the figure handle and the corrected coordinates and
% will be returned once the user presses the spacebar.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
%
% Andreas Horn

set(mcfig,'color','k');
axis off
setappdata(mcfig,'c_lims',[1800,2800]);
setappdata(mcfig,'selectrode',0);
setappdata(mcfig,'planecset',0);

set(mcfig,'KeyPressFcn',@ea_keystr);
set(mcfig, 'BusyAction','cancel', 'Interruptible','off');




setappdata(mcfig,'patientname',patientname);

%setappdata(mcfig,'trajectory',trajectory);

setappdata(mcfig,'origoptions',options); % store original options for further processing.


options.native=1;
setappdata(mcfig,'options',options);

if ~exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file')
    close(mcfig);
    msgbox('Please run pre-Reconstruct module first.');
    return
end

[coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

setappdata(mcfig,'origtrajectory',trajectory);
setappdata(mcfig,'manually_corrected',manually_corrected);

% initialize scene
updatescene([],[],mcfig);




% initialize toolbar
ht=uitoolbar(mcfig);
captions=getappdata(mcfig,'captions');
c_step=2;
minuscontrast=uipushtool(ht,'CData',ea_get_icn('contrastminus'),'TooltipString','Decrease Contrast [C]','ClickedCallback',{@setcontrast,'c',nan,mcfig});
pluscontrast=uipushtool(ht,'CData',ea_get_icn('contrastplus'),'TooltipString','Increase Contrast [V]','ClickedCallback',{@setcontrast,'v',nan,mcfig});
minusoffset=uipushtool(ht,'CData',ea_get_icn('extleft'),'TooltipString','Increase Offset [B]','ClickedCallback',{@setcontrast,'b',nan,mcfig});
plusoffset=uipushtool(ht,'CData',ea_get_icn('extright'),'TooltipString','Decrease Offset [N]','ClickedCallback',{@setcontrast,'n',nan,mcfig});

eltog(1)=uitoggletool(ht,'CData',ea_get_icn('el0'),'TooltipString','Select Electrode 0 [0]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(2)=uitoggletool(ht,'CData',ea_get_icn('el3'),'TooltipString','Select Electrode 3 [3]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(3)=uitoggletool(ht,'CData',ea_get_icn('el4'),'TooltipString','Select Electrode 4 [4]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(4)=uitoggletool(ht,'CData',ea_get_icn('el7'),'TooltipString','Select Electrode 7 [7]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});

rightview=uipushtool(ht,'CData',ea_get_icn('elR'),'TooltipString','Set view from Right [R]','ClickedCallback',{@ea_view,'r'});
leftview=uipushtool(ht,'CData',ea_get_icn('elL'),'TooltipString','Set view from Left [L]','ClickedCallback',{@ea_view,'l'});
antview=uipushtool(ht,'CData',ea_get_icn('elA'),'TooltipString','Set view from Anterior [A]','ClickedCallback',{@ea_view,'a'});
postview=uipushtool(ht,'CData',ea_get_icn('elP'),'TooltipString','Set view from Posterior [P]','ClickedCallback',{@ea_view,'p'});

%xview=uipushtool(ht,'CData',ea_get_icn('elX'),'TooltipString','Set view from X-Direction [X]','ClickedCallback',{@ea_view,'x'});
%yview=uipushtool(ht,'CData',ea_get_icn('elY'),'TooltipString','Set view from Y-Direction [Y]','ClickedCallback',{@ea_view,'y'});

% rotleft=uipushtool(ht,'CData',ea_get_icn('rotleft'),'TooltipString','Rotate Electrode counter-clockwise','ClickedCallback',{@ea_rotate,'cc',mcfig});
% rotright=uipushtool(ht,'CData',ea_get_icn('rotright'),'TooltipString','Rotate Electrode clockwise','ClickedCallback',{@ea_rotate,'c',mcfig});

rotleftcw=uipushtool(ht,'CData',ea_get_icn('rotleftcw'),'TooltipString','Rotate left electrode clockwise','ClickedCallback',{@ea_rotate,'lc',mcfig});
rotleftccw=uipushtool(ht,'CData',ea_get_icn('rotleftccw'),'TooltipString','Rotate left electrode counterclockwise','ClickedCallback',{@ea_rotate,'lcc',mcfig});
rotrightcw=uipushtool(ht,'CData',ea_get_icn('rotrightcw'),'TooltipString','Rotate right electrode clockwise','ClickedCallback',{@ea_rotate,'rc',mcfig});
rotrightccw=uipushtool(ht,'CData',ea_get_icn('rotrightccw'),'TooltipString','Rotate right electrode counterclockwise','ClickedCallback',{@ea_rotate,'rcc',mcfig});

mni=uitoggletool(ht,'CData',ea_get_icn('mninative'),'TooltipString','Toggle MNI vs. Native space','State','off','OnCallback',{@updatescene,mcfig,'mni'},'OffCallback',{@updatescene,mcfig,'native'});

finish_mc=uipushtool(ht,'CData',ea_get_icn('done'),'TooltipString','Finish manual corrections [space]','ClickedCallback',{@robotSpace});

captoggle=uitoggletool(ht,'CData',ea_get_icn('labels'),'TooltipString','Orientation','OnCallback',{@objvisible,captions},'OffCallback',{@objinvisible,captions},'State','on');



setappdata(mcfig,'eltog',eltog);


elplot=getappdata(mcfig,'elplot');

try
    realcoords_plot=getappdata(mcfig,'realcoords_plot');
end
trajectory_plot=getappdata(mcfig,'trajectory_plot');
planes=getappdata(mcfig,'planes');



%% Manual height correction here:
%set(mcfig,'Position',[10 400 700 500])

drawnow;

disp('Manual correction: Press arrows to adjust, space to end adjustment. For more shortcuts hover over buttons in the menu-bar.');


%% export variables to figure
setappdata(mcfig,'eltog',eltog);
setappdata(mcfig,'markers',markers);



function ea_endfcn(mcfig)
% This subfunction terminates the process of manual correction and saves
% results.
    disp('Saving results.');

    ea_busyaction('on',mcfig,'reco');
%markers=getappdata(gcf,'markers');
%trajectory=getappdata(gcf,'trajectory');
    options=getappdata(mcfig,'origoptions');

    options.hybridsave=1;
    options.native=1;
    [coords_mm,trajectory,markers,elmodel]=ea_load_reconstruction(options);

    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);
    options=getappdata(mcfig,'origoptions');
    try options=rmfield(options,'hybridsave'); end
ea_busyaction('off',mcfig,'reco');
close(mcfig)


% save results

%coords_mm=ea_resolvecoords(markers,options);
%elmodel=options.elmodel;

%    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);

%disp('Done.');

if options.autoimprove
    disp('Storing results in template.');
    ea_export_templates(coords_mm{1}(1:4,:),trajectory{1},options.patientname,options,'r')
    ea_export_templates(coords_mm{2}(1:4,:),trajectory{2},options.patientname,options,'l')
    disp('Done.');
end

%% methods dump:
ea_methods(options,...
            ['DBS-Electrodes were manually localized based on post-operative acquisitions using a tool specifically designed for this task (as implemented in Lead-DBS software',...
            '; Horn & Kuehn 2005; SCR_002915; http://www.lead-dbs.org).'],...
            {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});


% continue with rest of the program schedule..

ea_write(options);

% Callback invoked when user presses the 'Update all' button in the 'Manual electrodes head / tail setting' figure.
function manualMarkerSettingUpdateAllCb(source,event)
    manualMarkerSettingCbImpl(source,'all');
% Callback invoked when user presses the 'Cancel' button in the 'Manual electrodes head / tail setting' figure.
function manualMarkerSettingCancelCb(source,event)
    closereq();
% Callbacks invoked when user presses 'Update R0/R3/L3/L0' buttons in the 'Manual electrodes head / tail setting' figure.
function manualMarkerSettingUpdateR0Cb(source,event)
    manualMarkerSettingCbImpl(source,'R0');
function manualMarkerSettingUpdateR3Cb(source,event)
    manualMarkerSettingCbImpl(source,'R3');
function manualMarkerSettingUpdateL0Cb(source,event)
    manualMarkerSettingCbImpl(source,'L0');
function manualMarkerSettingUpdateL3Cb(source,event)
    manualMarkerSettingCbImpl(source,'L3');

function manualMarkerSettingCbImpl(source,action)
    userData=get(get(source,'parent'),'UserData');
    mcfig=userData.mcfig;
    options=getappdata(mcfig,'options');
    elplot=getappdata(mcfig,'elplot');
    mplot=getappdata(mcfig,'mplot');

    % load the current parameters
    [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
    % update markers according to manually specified positions
    mRH=str2num(get(userData.tvRH,'String'));
    mRT=str2num(get(userData.tvRT,'String'));
    mLH=str2num(get(userData.tvLH,'String'));
    mLT=str2num(get(userData.tvLT,'String'));
    switch action
        case 'all'
            markers(1).head=mRH;
            markers(1).tail=mRT;
            markers(2).head=mLH;
            markers(2).tail=mLT;
        case 'R0'
            markers(1).tail=markers(1).tail+(mRH-markers(1).head);
            markers(1).head=mRH;
        case 'R3'
            markers(1).head=markers(1).head+(mRT-markers(1).tail);
            markers(1).tail=mRT;
        case 'L0'
            markers(2).tail=markers(2).tail+(mLH-markers(2).head);
            markers(2).head=mLH;
        case 'L3'
            markers(2).head=markers(2).head+(mLT-markers(2).tail);
            markers(2).tail=mLT;
        otherwise
            error(['unknown action "' action '"']);
    end

    % close the 'Manual electrodes head / tail setting' figure
    closereq();

    % save the parameters
    if isfield(options,'hybridsave')
        options=rmfield(options,'hybridsave');
    end
    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);

    % the rest will be done by the code that created the 'Manual electrodes head / tail setting'
    % figure in the 'ea_keystr' function

function ea_keystr(mcfig,event)
%    pause
%commnd=get (gcf, 'CurrentKey');


%% get vars
eltog=getappdata(mcfig,'eltog');
elplot=getappdata(mcfig,'elplot');
mplot=getappdata(mcfig,'mplot');
options=getappdata(mcfig,'options');

if ~exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file')
    close(mcfig);
    return
end
[coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);



commnd=event.Character;
switch lower(commnd)
    case ' '
        %markers=getappdata(mcfig,'markers');
        ea_endfcn(mcfig);
        return
    case {'x','a','p','y','l','r'} % view angles.
        %markers=getappdata(mcfig,'markers');
        ea_view(nan,nan,commnd);

    case {'0','3','4','7'}
        switch lower(commnd)
            case '0'
                selectrode=1;
            case '3'
                selectrode=2;
            case '4'
                selectrode=3;
            case '7'
                selectrode=4;
        end
        oselectrode=getappdata(mcfig,'selectrode');
        if selectrode==oselectrode % toggle had already been clicked -> deselect all.
            % reset all toggletools
            for i=1:4
                set(eltog(i),'State','off');
            end
            setappdata(mcfig,'selectrode',0);
            updatescene([],[],mcfig);
        else
            % clear all toggletools.
            for i=1:4
                set(eltog(i),'State','off');
            end

            % set the correct toggletool again.
            set(eltog(selectrode),'State','on');

            % store selected electrode in appdata.
            setappdata(mcfig,'selectrode',selectrode);
            updatescene([],[],mcfig);
        end
        clear oselectrode
    case {'c','v','b','n'}
        setcontrast(nan,nan,event.Character,event.Modifier,mcfig);
%     case 'v'
%         increasecontrast(nan,nan,event);
%     case 'b'
%         increaseoffset(nan,nan,event);
%     case 'n'
%         decreaseoffset(nan,nan,event);

    case {'s','d'}

    case 'm' % manual electrodes head / tail setting
        % get screen size
        scrSize=get(0,'ScreenSize');
        % create a small figure in the bottom part of the screen
        FIG_WIDTH=.4;
        f=figure('Name','Manual electrodes head / tail setting','Numbertitle','off',...
            'MenuBar','none',...
            'Position',[scrSize(1)+scrSize(3)*(.5-FIG_WIDTH/2) scrSize(2) scrSize(3)*FIG_WIDTH scrSize(4)/5]);
        B=.05; % elementary border
        H=.18-B; % height of each element
        % TI = text infos
        XTI=B;
        WTI=.15; % width
        % TV = (editable) text values
        WTV=.35;
        XTV=XTI+WTI+B;
        % BT = buttons on the right side
        WBT=1-WTI-WTV-4*B;
        XBT=XTV+WTV+B;
        %
        YRH=1-.2-B;
        YRT=YRH-H-B;
        YLH=YRT-H-B;
        YLT=YLH-H-B;
        YBT_ALL=YLT-H-B;
        if ismember(1,options.sides)
            % right side
            tiRH=uicontrol('Style','text','String','R0 (Head):','HorizontalAlignment','right',...
                'Units','normalized','Position',[XTI YRH WTI H]);
            tvRH=uicontrol('Style','edit','String',num2str(markers(1).head,'%.4f '),...
                'HorizontalAlignment','center',...
                'Units','normalized','Position',[XTV YRH WTV H]);
            moveR0=uicontrol('Style','pushbutton','String','Update R0, Move R3 accordingly',...
                'Units','normalized','Position',[XBT YRH WBT H],'Callback',@manualMarkerSettingUpdateR0Cb);
            tiRT=uicontrol('Style','text','String','R3 (Tail):','HorizontalAlignment','right',...
                'Units','normalized','Position',[XTI YRT WTI H]);
            tvRT=uicontrol('Style','edit','String',num2str(markers(1).tail,'%.4f '),...
                'HorizontalAlignment','center',...
                'Units','normalized','Position',[XTV YRT WTV H]);
            moveR3=uicontrol('Style','pushbutton','String','Update R3, Move R0 accordingly',...
                'Units','normalized','Position',[XBT YRT WBT H],'Callback',@manualMarkerSettingUpdateR3Cb);
        end
        if ismember(2,options.sides)
            % left side
            tiLH=uicontrol('Style','text','String','L0 (Head):','HorizontalAlignment','right',...
                'Units','normalized','Position',[XTI YLH WTI H]);
            tvLH=uicontrol('Style','edit','String',num2str(markers(2).head,'%.4f '),...
                'HorizontalAlignment','center',...
                'Units','normalized','Position',[XTV YLH WTV H]);
            moveL0=uicontrol('Style','pushbutton','String','Update L0, Move L3 accordingly',...
                'Units','normalized','Position',[XBT YLH WBT H],'Callback',@manualMarkerSettingUpdateL0Cb);
            tiLT=uicontrol('Style','text','String','L3 (Tail):','HorizontalAlignment','right',...
                'Units','normalized','Position',[XTI YLT WTI H]);
            tvLT=uicontrol('Style','edit','String',num2str(markers(2).tail,'%.4f '),...
                'HorizontalAlignment','center',...
                'Units','normalized','Position',[XTV YLT WTV H]);
            moveL3=uicontrol('Style','pushbutton','String','Update L3, Move L0 accordingly',...
                'Units','normalized','Position',[XBT YLT WBT H],'Callback',@manualMarkerSettingUpdateL3Cb);
        end
        updateButton=uicontrol('Style','pushbutton','String','Update all',...
            'Units','normalized','Position',[XTV+B YBT_ALL WTV-2*B H],'Callback',@manualMarkerSettingUpdateAllCb);
        cancelButton=uicontrol('Style','pushbutton','String','Cancel',...
            'Units','normalized','Position',[XBT+B YBT_ALL WBT-2*B H],'Callback',@manualMarkerSettingCancelCb);
        % store the UI elements in user data property of the new figure
        userData=struct('tvRH',tvRH,'tvRT',tvRT,'tvLH',tvLH,'tvLT',tvLT,'mcfig',mcfig);
        set(f,'UserData',userData);

        % wait for the figure closing (either by pressing the bu button, or by closing it)
        waitfor(f)

        % update the scene
        updatescene([],[],mcfig);

        % reload the parameters
        [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

    otherwise % arrow keys, plus, minus

        if ismember(event.Key,{'rightarrow','leftarrow','uparrow','downarrow'}) || ismember(event.Character,{'+','-','*','_'})
        selectrode=getappdata(mcfig,'selectrode');
        if ~selectrode % no electrode is highlighted, move electrodes alongside trajectory or increase/decrease spacing.
            %markers=getappdata(mcfig,'markers');
            %trajectory=getappdata(mcfig,'trajectory');
            [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);


            markers=ea_correctcoords(markers,trajectory,event);
            if isfield(options,'hybridsave')
                options=rmfield(options,'hybridsave');
            end
            ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);

            %            setappdata(mcfig,'markers',markers);
            updatescene([],[],mcfig);
            [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

            %           markers=getappdata(mcfig,'markers');
        else % electrode is highlighted. Move in xy dirs.

%             markers=getappdata(mcfig,'markers');
%             trajectory=getappdata(mcfig,'trajectory');
            [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
            movedcoords=moveonecoord(markers,selectrode,event); % move the correct coord to the correct direction.

            for side=options.sides
                    set(mplot(1,side),'XData',movedcoords(side).head(1),'YData',movedcoords(side).head(2),'ZData',movedcoords(side).head(3))
                    set(mplot(2,side),'XData',movedcoords(side).tail(1),'YData',movedcoords(side).tail(2),'ZData',movedcoords(side).tail(3))
            end
%            setappdata(mcfig,'markers',markers);
            ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);

            updatescene([],[],mcfig);
            %markers=getappdata(mcfig,'markers');
            [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

            %update_coords(elplot(selectrode),markers,trajectory,movedcoords); % refresh scene view (including update for all other electrodes).
            %setappdata(gcf,'markers',markers);
            %setappdata(mcfig,'trajectory',trajectory);

        end
        end
end
%end
cnt=1;
options=getappdata(mcfig,'options');
coords_mm=ea_resolvecoords(markers,options);
           % [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

for side=options.sides
    set(mplot(1,side),'XData',markers(side).head(1),'YData',markers(side).head(2),'ZData',markers(side).head(3));
    set(mplot(2,side),'XData',markers(side).tail(1),'YData',markers(side).tail(2),'ZData',markers(side).tail(3));
    for el=1:length(elplot)/2
        set(elplot(cnt),'XData',coords_mm{side}(el,1),'YData',coords_mm{side}(el,2),'ZData',coords_mm{side}(el,3))
        cnt=cnt+1;
    end
end
refreshdata(elplot,'caller')
drawnow










function hdtrajectory=genhd_inside(trajectory)

resolution=20;

hdtrajectory(:,1)=interp1q([1:length(trajectory)]',trajectory(:,1),[1:1/resolution:length(trajectory)]');
hdtrajectory(:,2)=interp1q([1:length(trajectory)]',trajectory(:,2),[1:1/resolution:length(trajectory)]');
hdtrajectory(:,3)=interp1q([1:length(trajectory)]',trajectory(:,3),[1:1/resolution:length(trajectory)]');



function updatescene(varargin)

hobj=varargin{1};
ev=varargin{2};
mcfig=varargin{3};
firstrun=getappdata(mcfig,'firstrun');
contrast=getappdata(mcfig,'contrast');
offset=getappdata(mcfig,'offset');
if isempty(contrast), contrast=0.8; end
if isempty(offset), offset=0.3; end
setappdata(mcfig,'offset',offset); setappdata(mcfig,'contrast',contrast);
%% inputs:
options=getappdata(mcfig,'options');

patientname=getappdata(mcfig,'patientname');
%markers=getappdata(mcfig,'markers');

if nargin==4
        ea_busyaction('on',gcf,'reco');

    options.hybridsave=1;
    [coords_mm,trajectory,markers,elmodel]=ea_load_reconstruction(options);
    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);
    options=rmfield(options,'hybridsave');
    space=varargin{4};
        ea_busyaction('del',gcf,'reco');

else
    space=getappdata(mcfig,'space');
    if isempty(space)
        space='native';
    end
end
setappdata(mcfig,'space',space);
switch space
    case 'mni'
        options.native=0;
    case 'native'
        options.native=1;
end
setappdata(mcfig,'options',options);


[~,~,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);



if isempty(firstrun) && ~manually_corrected % resize electrode to default spacing.
    [~,trajectory,markers]=ea_resolvecoords(markers,options,1);
    setappdata(mcfig,'firstrun',0);
else
    [~,trajectory,markers]=ea_resolvecoords(markers,options,0);
end

% for now, rotation will always be constant. This will be the place to
% insert rotation functions..

% for side=options.sides
%     normtrajvector=(markers(side).tail-markers(side).head)./norm(markers(side).tail-markers(side).head);
%     orth=null(normtrajvector)*(options.elspec.lead_diameter/2);
%     markers(side).x=markers(side).head+orth(:,1)';
%     markers(side).y=markers(side).head+orth(:,2)'; % corresponding points in reality
% end
for side=options.sides
    rotation=getappdata(gcf,'rotation');
    if isempty(rotation)
        rotation{1} = 0;
        rotation{2} = 0;
    end
    normtrajvector=(markers(side).tail-markers(side).head)./norm(markers(side).tail-markers(side).head);
    
    y(1) = -cos(0) * sin(deg2rad(rotation{side}));
    y(2) = (cos(0) * cos(deg2rad(rotation{side}))) + (sin(0) * sin(deg2rad(rotation{side})) * sin(0));
    y(3) = (-sin(0) * cos(deg2rad(rotation{side}))) + (cos(0) * sin(deg2rad(rotation{side})) * sin(0));
    y = y - (dot(y,normtrajvector) / (norm(normtrajvector) ^2)) * normtrajvector;
    x = cross(y,normtrajvector);
    
    markers(side).x=markers(side).head+x;
    markers(side).y=markers(side).head+y; % corresponding points in reality
end



%            [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);



%trajectory=getappdata(mcfig,'trajectory');
options=getappdata(mcfig,'options');
movedel=getappdata(mcfig,'movedel');
trajectory_plot=getappdata(mcfig,'trajectory_plot');
spacetext=getappdata(mcfig,'spacetext');
planes=getappdata(mcfig,'planes');
c_lims=getappdata(mcfig,'c_lims');

elplot = getappdata(mcfig,'elplot');
mplot = getappdata(mcfig,'mplot');

selectrode=getappdata(mcfig,'selectrode');

if ~isempty(selectrode) && selectrode>0
    coordhandle=mplot(selectrode);
end




movedmarkers=nan(4,3);
try
    movedmarkers(1,:)=[get(mplot(1),'xdata'),...
        get(mplot(1),'ydata'),...
        get(mplot(1),'zdata')];
end
try
    movedmarkers(2,:)=[get(mplot(2),'xdata'),...
        get(mplot(2),'ydata'),...
        get(mplot(2),'zdata')];
end
try
    movedmarkers(3,:)=[get(mplot(3),'xdata'),...
        get(mplot(3),'ydata'),...
        get(mplot(3),'zdata')];
end
try
    movedmarkers(4,:)=[get(mplot(4),'xdata'),...
        get(mplot(4),'ydata'),...
        get(mplot(4),'zdata')];
end

% xdata = cell2mat(get(mplot,'xdata'));
% ydata = cell2mat(get(mplot,'ydata'));
% zdata = cell2mat(get(mplot,'zdata'));
% movedmarkers=[xdata,ydata,zdata];

if selectrode
    
    [markers]=update_coords(coordhandle,markers,trajectory,movedmarkers,options);
end


[coords_mm,trajectory]=ea_resolvecoords(markers,options);




%% plot main figure

delete(trajectory_plot);
delete(spacetext);



mainax=subplot(4,5,[2:4,7:9,12:14,17:19]); % main plot
set(gca, 'LooseInset', [0,0,0,0]);
init=getappdata(mcfig,'init');
if isempty(init)
    view(0,0);
    axis off
    setappdata(mcfig,'init',1)
end


% delete prior captions
spacetext=getappdata(mcfig,'spacetext');
delete(spacetext);
captions=getappdata(mcfig,'captions');
delete(captions);

% Plot spacing distance info text and correct inhomogeneous spacings.
%emp_eldist(1)=mean([ea_pdist([markers(1).head;markers(1).tail]),ea_pdist([markers(2).head;markers(2).tail])])/3;
clear emp_eldist
for side=options.sides
    A{side}=sqrt(ea_sqdist(coords_mm{side}',coords_mm{side}'));
    emp_eldist{side}=sum(sum(tril(triu(A{side},1),1)))/(options.elspec.numel-1);
end
memp_eldist=mean([emp_eldist{:}]);

    [~,trajectory,markers]=ea_resolvecoords(markers,options,1,memp_eldist);


%% plot coords
hold on



if isempty(elplot) % first time plot electrode contacts
cnt=1;

    for side=options.sides
        mplot(1,side)=plot3(markers(side).head(1),markers(side).head(2),markers(side).head(3),'*','MarkerEdgeColor',[0.9 0.2 0.2],'MarkerFaceColor','none','MarkerSize',25);
        mplot(2,side)=plot3(markers(side).tail(1),markers(side).tail(2),markers(side).tail(3),'*','MarkerEdgeColor',[0.2 0.9 0.2],'MarkerFaceColor','none','MarkerSize',25);
        for el=1:size(coords_mm{side},1)
            elplot(cnt)=plot3(coords_mm{side}(el,1),coords_mm{side}(el,2),coords_mm{side}(el,3),'O','MarkerEdgeColor',[0.9 0.9 0.9],'MarkerFaceColor','none','MarkerSize',25);
            cnt=cnt+1;
        end
    end
    
   setappdata(mcfig,'elplot',elplot);
   setappdata(mcfig,'mplot',mplot);

else % update coordinates in elplot & mplot:
    cnt=1;

    for side=options.sides

        set(mplot(1,side),'XData',markers(side).head(1));
        set(mplot(1,side),'YData',markers(side).head(2));
        set(mplot(1,side),'ZData',markers(side).head(3));

        set(mplot(2,side),'XData',markers(side).tail(1));
        set(mplot(2,side),'YData',markers(side).tail(2));
        set(mplot(2,side),'ZData',markers(side).tail(3));

        for el=1:size(coords_mm{side},1)
            set(elplot(cnt),'XData',coords_mm{side}(el,1));
            set(elplot(cnt),'YData',coords_mm{side}(el,2));
            set(elplot(cnt),'ZData',coords_mm{side}(el,3));
            cnt=cnt+1;
        end
    end
    setappdata(mcfig,'elplot',elplot);
    setappdata(mcfig,'mplot',mplot);
end



try
    midpt=mean([markers(1).head;markers(2).head],1);
catch
    midpt=[0 0 0];
end


spacetext=text(midpt(1),midpt(2),midpt(3)-1,sprintf(['Electrode Spacing: ',num2str(memp_eldist),' mm\nRight Rotation: ',num2str(rotation{1}),' °\nLeft Rotation: ',num2str(rotation{2}),' °']),'Color','w','BackgroundColor','k','HorizontalAlignment','center');
set(mcfig,'name',[options.patientname,', Electrode Spacing: ',num2str(memp_eldist),' mm.']);
setappdata(mcfig,'spacetext',spacetext);

% %% plot trajectory lines

for side=1:length(options.sides)

    try

        if ~isempty(trajectory{side})

            if options.verbose>1; 
                trajectory_plot(side)=plot3(trajectory{side}(:,1),trajectory{side}(:,2),trajectory{side}(:,3),'color',[0.3,0.5,0.9],'linew',1.5); 
            end

        end
    end
end


delete(planes);
clear planes
planecnt=1;
%% plot slices in x and y planes

for doxx=0:1
    for side=options.sides
        %try
            sample_width=20-doxx*5; % a bit smaller sample size in x direction to avoid overlap.
            meantrajectory=genhd_inside(trajectory{side});
            clear imat
            %% sample plane left and right from meantrajectory

            if doxx
                Vcor=getV(mcfig,'Vcor',options);
                imat=ea_resample_planes(Vcor,meantrajectory',sample_width,doxx,0.2);

            else
                Vsag=getV(mcfig,'Vsag',options);
                imat=ea_resample_planes(Vsag,meantrajectory',sample_width,doxx,0.2);

            end


            colormap gray


            if doxx % span surface in x direction
                spanvector=[sample_width,0,0];
            else % span surface in y direction
                spanvector=[0,sample_width,0];
            end

            boundingbox=[meantrajectory(1,:)-spanvector;...
                meantrajectory(1,:)+spanvector;...
                meantrajectory(end,:)-spanvector;...
                meantrajectory(end,:)+spanvector];


            xx=[boundingbox(1,1),boundingbox(2,1);boundingbox(3,1),boundingbox(4,1)];
            yy=[boundingbox(1,2),boundingbox(2,2);boundingbox(3,2),boundingbox(4,2)];
            zz=[boundingbox(1,3),boundingbox(2,3);boundingbox(3,3),boundingbox(4,3)];


            alphamap=imat;
            alphamap(:)=0.9;

            if ~getappdata(mcfig,'planecset') % initially and once set contrast based on image data.

                if options.modality==1 % MR
                    c_lims=[ea_nanmean(imat(:))-ea_nanstd(imat(:))-3*ea_nanstd(imat(:)),ea_nanmean(imat(:))-ea_nanstd(imat(:))+3*ea_nanstd(imat(:))];
                elseif options.modality==2 % CT

                    lthresh=800; % initial guesses for CT
                    uthresh=2800;
                    try % try estimating a better guess..
                        for tries=1:5
                            timat=imat;
                            timat(timat<lthresh)=0;
                            timat(timat>uthresh)=0;

                            nomi=ea_nmi(round(imat),round(timat));
                            if nomi>0.9
                                break
                            else
                                lthresh=lthresh+randn(1)*200;
                                uthresh=uthresh+randn(1)*200;
                                if lthresh>=uthresh
                                    lthresh=uthresh-500;
                                end
                            end
                        end
                    end
                   % disp(['Lthresh: ',num2str(lthresh),'; Uthresh: ',num2str(uthresh),'.']);
                    c_lims=[lthresh,uthresh]; % Initial guess, CT


                end
                %caxis(c_lims);
                caxis([0,1]);
                setappdata(mcfig,'c_lims',c_lims);
                setappdata(mcfig,'planecset',1);
            end

imat=ea_contrast(imat,contrast,offset);

            planes(planecnt)=surface('XData',xx,'YData',yy,'ZData',zz,'CData',imat,'alphadata',alphamap,'FaceAlpha', 'texturemap','FaceColor','texturemap','EdgeColor','none','alphadatamapping','none');

            planecnt=planecnt+1;
            if ~doxx

                if side==1

%                     captions(1)=text(0,... % x
%                         ((min(boundingbox(:,2))+max(boundingbox(:,2)))/2)+20,... % y
%                         0,... % z
%                         'A','Color','w','BackgroundColor','k');
%                     captions(2)=text(0,... % x
%                         ((min(boundingbox(:,2))+max(boundingbox(:,2)))/2)-20,... % y
%                         0,... % z
%                         'P','Color','w','BackgroundColor','k');
%                     % captions(1)=text(0,0,0,... % z
%                     %  'C','Color','w');
%
%                     captions(3)=text(40,... % x
%                         0,... % y
%                         0,... % z
%                         'R','Color','w','BackgroundColor','k');
%
%                     captions(4)=text(-40,... % x
%                         0,... % y
%                         0,... % z
%                         'L','Color','w','BackgroundColor','k');
%


captions(1)=text(midpt(1),... % x
    midpt(2)+10,... % y
    midpt(3)+10,... % z
    'A','Color','w','BackgroundColor','k');
captions(2)=text(midpt(1),... % x
    midpt(2)-10,... % y
    midpt(3)+10,... % z
    'P','Color','w','BackgroundColor','k');
captions(3)=text(midpt(1)-10,... % x
    midpt(2),... % y
    midpt(3)+10,... % z
    'L','Color','w','BackgroundColor','k');
captions(4)=text(midpt(1)+10,... % x
    midpt(2),... % y
    midpt(3)+10,... % z
    'R','Color','w','BackgroundColor','k');



                    setappdata(mcfig,'captions',captions);

                end
            else

            end
        %end
    end
end
%caxis([c_lims(1) c_lims(2)]);
caxis([0,1]);



%% plot axial planes on the right hand side of the figure


Vtra=getV(mcfig,'Vtra',options);

mks=nan(4,3); % always assign 4 markers, no matter if only right or left electrode selected. fill with nans
try mks(1,:)=markers(1).head; end
try mks(2,:)=markers(1).tail; end
try mks(3,:)=markers(2).head; end
try mks(4,:)=markers(2).tail; end

        mks=Vtra.mat\[mks,ones(size(mks,1),1)]';
        mks=mks(1:3,:)';


%title(['Electrode ',num2str(el-1),', transversal view.']);
wsize=10;
cmap=[1,4,5,8];

for subpl=getsuplots(options.sides)
    
    subplot(4,5,subpl*5)

    slice=ea_sample_slice(Vtra,'tra',wsize,'vox',mks,subpl);
    slice=ea_contrast(slice,contrast,offset);
    try
        imagesc(slice,[ea_nanmean(slice(slice>0))-3*ea_nanstd(slice(slice>0)) ea_nanmean(slice(slice>0))+3*ea_nanstd(slice(slice>0))]);
    catch
        imagesc(slice);
    end

    hold on



    if selectrode && subpl==selectrode

        fc='y';
    else
        if ismember(subpl,[1,3])
        fc='r';
        else
            fc='g';
        end
    end

    warnStruct = warning('off','MATLAB:hg:willberemoved');
    plot((wsize+1)*2,(wsize+1)*2,'*','MarkerSize',15,'MarkerEdgeColor',fc,'LineWidth',2,'LineSmoothing','on');
    warning(warnStruct);
    hold off
    axis square
    axis off
    %caxis([c_lims(1) c_lims(2)]);
caxis([0,1]);
end



%% plot electrode model to the left side (static)
legplot=getappdata(mcfig,'legplot');
if isempty(legplot)
elax=subplot(4,5,[1,6,11,16]); % left electrode plot
axis off
load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname])

% visualize
cnt=1;
X=eye(4);
hold on
for ins=1:length(electrode.insulation)
    electrode.insulation(ins).vertices=X*[electrode.insulation(ins).vertices,ones(size(electrode.insulation(ins).vertices,1),1)]';
    electrode.insulation(ins).vertices=electrode.insulation(ins).vertices(1:3,:)';
    elrender{side}(cnt)=patch(electrode.insulation(ins));
    ea_specsurf(elrender{side}(cnt),options.elspec.lead_color,0.5);
    cnt=cnt+1;
end
for con=1:length(electrode.contacts)
    electrode.contacts(con).vertices=X*[electrode.contacts(con).vertices,ones(size(electrode.contacts(con).vertices,1),1)]';
    electrode.contacts(con).vertices=electrode.contacts(con).vertices(1:3,:)';
    elrender{side}(cnt)=patch(electrode.contacts(con));

    ea_specsurf(elrender{side}(cnt),options.elspec.contact_color,0.5);

    cnt=cnt+1;
end

plot3(electrode.head_position(1),electrode.head_position(2),electrode.head_position(3),'*r','MarkerSize',15)
plot3(electrode.tail_position(1),electrode.tail_position(2),electrode.tail_position(3),'*g','MarkerSize',15)
axis([-2,2,-2,2,0,16])
set(elax,'XLimMode','manual'),set(elax,'YLimMode','manual'),set(elax,'ZLimMode','manual')
axis manual
axis equal
view(0,0);

%light('Position',[0 -5 10]);

text(0,0,14,options.elmodel,'color','w');


setappdata(mcfig,'legplot',1);

end


%% outputs

% try
%     setappdata(resultfig,'realcoords_plot',realcoords_plot);
% end

set(mcfig,'CurrentAxes',mainax);
setappdata(mcfig,'planes',planes);







%% outputs

setappdata(mcfig,'elplot',elplot);
setappdata(mcfig,'mplot',mplot);
setappdata(mcfig,'movedel',movedel);

if isfield(options,'hybridsave')
    options=rmfield(options,'hybridsave');
end

    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);


%setappdata(mcfig,'markers',markers);
% try
%     setappdata(resultfig,'realcoords_plot',realcoords_plot);
% end
setappdata(mcfig,'trajectory_plot',trajectory_plot);
setappdata(mcfig,'planes',planes);
%setappdata(mcfig,'trajectory',trajectory);

function sp=getsuplots(sides)
if isequal(sides,[1:2])
    sp=1:4;
elseif isequal(sides,1)
    sp=1:2;
elseif isequal(sides,2)
    sp=3:4;
end

function V=getV(mcfig,ID,options)
if options.native
    addon='_unnormalized';
else
    addon='';
end

switch options.modality
    case 1 % MR
V=getappdata(mcfig,[ID,addon]);

if isempty(V)
%     flags.interp=4;
%     flags.wrap=[0,0,0];
%         d       = [flags.interp*[1 1 1]' flags.wrap(:)];
     switch ID
        case 'Vcor'
            try
                V=spm_vol([options.root,options.patientname,filesep,options.prefs.(['cornii',addon])]);
            catch
                V=spm_vol([options.root,options.patientname,filesep,options.prefs.(['tranii',addon])]);
            end
        case 'Vtra'

            V=spm_vol([options.root,options.patientname,filesep,options.prefs.(['tranii',addon])]);

        case 'Vsag'
            try
                V=spm_vol([options.root,options.patientname,filesep,options.prefs.(['sagnii',addon])]);


            catch
                try
                    V=spm_vol([options.root,options.patientname,filesep,options.prefs.(['cornii',addon])]);
                catch
                    V=spm_vol([options.root,options.patientname,filesep,options.prefs.(['tranii',addon])]);
                end
            end
    end

 %   C=spm_bsplinc(V,d);
end
setappdata(mcfig,[ID,addon],V);

    case 2 % CT - ignore wishes, always feed out V as CT.
        if options.native
            V=getappdata(mcfig,'VCTnative');
            if isempty(V)
                V=spm_vol([options.root,options.patientname,filesep,options.prefs.ctnii_coregistered]);
                setappdata(mcfig,'VCTnative',V);
            end
        else
            V=getappdata(mcfig,'VCTmni');
            if isempty(V)
                V=spm_vol([options.root,options.patientname,filesep,options.prefs.ctnii]);
                setappdata(mcfig,'VCTmni',V);
            end
        end
end
%setappdata(mcfig,['C',ID,addon],C);



function movedel=whichelmoved(coordhandle)

mplot = getappdata(gcf,'mplot');
for side=1:2
    for el=1:2;

        if coordhandle==mplot(el,side)
            movedel=sub2ind([2,2],el,side);
        end

    end
end

function   [numarkers,nutrajectory]=update_coords(coordhandle,markers,trajectory,movedmarkers,options)

% movedcoords are NOT in cell format!

movel=whichelmoved(coordhandle);

mplot=getappdata(gcf,'mplot');
origtrajectory=getappdata(gcf,'origtrajectory');
%whichonemoved=1;

numarkers=markers;
nutrajectory=trajectory;

switch movel

    case 1 % lower right
        if ismember(1,options.sides)
            side=1;
        olddist=abs(ea_pdist([markers(1).head;markers(1).tail]));
        refel=4;
        changecoord=1:3;
        cchangecoord=1;
        usefit=1;
        spin=-1;

        refpt=markers(1).tail;
        movedpt=movedmarkers(1,:);
        move='head';
        end
    case 3 % lower left
        if ismember(2,options.sides)
            side=2;
        olddist=abs(ea_pdist([markers(2).head;markers(2).tail]));
        refel=4;
        changecoord=1:3;
        cchangecoord=2;
        usefit=2;
        spin=-1;

        refpt=markers(2).tail;
        try
        movedpt=movedmarkers(3,:);
        catch
            keyboard
        end
        move='head';
        end
    case 2 % upper right
        if ismember(1,options.sides)
            side=1;
            olddist=abs(ea_pdist([markers(1).head;markers(1).tail]));
            refel=1;
            changecoord=2:4;
            cchangecoord=1;
            usefit=1;
            spin=1;
            
            refpt=markers(1).head;
            movedpt=movedmarkers(2,:);
            move='tail';
        end
    case 4 % upper left
        if ismember(2,options.sides)
            side=2;
        olddist=abs(ea_pdist([markers(2).head;markers(2).tail]));
        refel=1;
        changecoord=2:4;
        cchangecoord=2;

        usefit=2;
        spin=1;

        refpt=markers(2).head;
        movedpt=movedmarkers(4,:);
         move='tail';
        end
end


helppt=[movedpt(1:2),refpt(3)]; % this point is constructed as a projection of the movedpoint to the z-axis at height of refpt.
zdistfromhelppt=sqrt(abs((olddist)^2-(ea_pdist([helppt;refpt]))^2)); % use Pythagorean theorem to calculate distance on z-axis.
% here, 3*olddist is the hypotenuse, ea_pdist between the helper point and the
% reference point is one leg of the triangle, zdistfromhelppt is the
% leg that is calculated (the distance from the helper point to the new
% moved point (fulfilling the prerequisite that new moved point is at the
% same distance than the old point).
movedpt=helppt+spin*[0 0 zdistfromhelppt]; % set movedpt to be equidistant.

trajpts=[refpt;movedpt];

nutraj=diff(trajpts);
nutraj=nutraj/norm(nutraj);
%disp(num2str(nutraj));

% generate new coords


try
xc=refpt+(nutraj*(olddist));
numarkers(cchangecoord)=setfield(numarkers(cchangecoord),move,xc);
    set(mplot(cchangecoord,side),'xdata',xc(1));
    set(mplot(cchangecoord,side),'ydata',xc(2));
    set(mplot(cchangecoord,side),'zdata',xc(3));
catch
    keyboard
end

% generate new trajectory
nutrajectory{usefit}=[];
cnt=1;

for ix=-50:0.5:50

    thispoint=refpt+(nutraj*ix);
    if thispoint(3)<max(origtrajectory{usefit}(:,3)) && thispoint(3)>min(origtrajectory{usefit}(:,3))
        nutrajectory{usefit}(cnt,:)=thispoint;
        cnt=cnt+1;
    end


end
nutrajectory{usefit}=sortrows(nutrajectory{usefit},-3); % Make sure that trajectory is always listed in correct order.





function markers=moveonecoord(markers,selectrode,command)

grone=[0.1,1]; % step-sizes

movedcoords=[markers(1).head;markers(1).tail;markers(2).head;markers(2).tail];

switch command.Key
    case 'leftarrow'
        movedcoords(selectrode,1)=movedcoords(selectrode,1)-grone(1+ismember('shift',command.Modifier));
    case 'rightarrow'
        movedcoords(selectrode,1)=movedcoords(selectrode,1)+grone(1+ismember('shift',command.Modifier));
    case 'uparrow'
        movedcoords(selectrode,2)=movedcoords(selectrode,2)+grone(1+ismember('shift',command.Modifier));

    case 'downarrow'
        movedcoords(selectrode,2)=movedcoords(selectrode,2)-grone(1+ismember('shift',command.Modifier));



end
markers(1).head=movedcoords(1,:);
markers(1).tail=movedcoords(2,:);
markers(2).head=movedcoords(3,:);
markers(2).tail=movedcoords(4,:);


function objvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);

function objinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);

function ea_rotate(hobj,ev,ccw,mcfig)
rotation=getappdata(gcf,'rotation'); % rotation angle in degrees
if isempty(rotation)
    rotation{1} = 0;
    rotation{2} = 0;
end
switch ccw
    case 'rc'
        rotation{1}=rotation{1}-15;
    case 'rcc'
        rotation{1}=rotation{1}+15;
    case 'lc'
        rotation{2}=rotation{2}-15;
    case 'lcc'
        rotation{2}=rotation{2}+15;        
end
setappdata(gcf,'rotation',rotation);
updatescene([],[],mcfig);

function setcontrast(hobj,ev,key,modifier,mcfig)
% c_lims=getappdata(gcf,'c_lims');
% comms={'v','c','b','n'}; % key commands
% perfs=[1 -1 % actions to perform on key commands
%       -1 1
%        1 1
%       -1 -1];
% kern=(c_lims(2)-c_lims(1))/20; % gain to correct contrast
% 
 doshift=any(ismember('shift',modifier));
 
% 
% c_lims=c_lims+(perfs(ismember(comms,lower(key)),:)*(kern*(doshift+1)));
% setappdata(gcf,'c_lims',c_lims);
doshift=1+doshift*4;
% new code:
contrast=getappdata(mcfig,'contrast');
offset=getappdata(mcfig,'offset');
grain=0.1;

switch lower(key)
    case 'v'
        contrast=contrast+grain*doshift;
    case 'c'
        contrast=contrast-grain*doshift;
    case 'b'
        offset=offset-grain*doshift*2;
    case 'n'
        offset=offset+grain*doshift*2;
end
% if offset<=0
%     offset=0;
% end
if contrast<0.1;
    contrast=0.1;
end
setappdata(mcfig,'contrast',contrast);
setappdata(mcfig,'offset',offset);
%disp(['Contrast: ',num2str(contrast),', Offset: ',num2str(offset),'.']);
updatescene([],[],mcfig);


function selectelectrode(hobj,ev)
% reset all toggletools

warning ('off','all');

eltog=getappdata(gcf,'eltog');
elplot=getappdata(gcf,'elplot');
mplot=getappdata(gcf,'mplot');
options=getappdata(gcf,'options');
    set(mplot(1,options.sides),'MarkerEdgeColor','r');
    set(mplot(2,options.sides),'MarkerEdgeColor','g');

for i=getsuplots(options.sides)
    set(eltog(i),'State','off');
    if eltog(i)==hobj;
        selectrode=i;
    end
end
if ~exist('selectrode','var')
    return
end

% set the clicked toggletool again.
set(eltog(selectrode),'State','on');

set(mplot(selectrode),'MarkerEdgeColor','y');

% store selected electrode in appdata.
setappdata(gcf,'selectrode',selectrode);
warning ('off','all');




function deselectelectrode(hobj,ev)
% reset all toggletools
mplot=getappdata(gcf,'mplot');
options=getappdata(gcf,'options');
    set(mplot(1,options.sides),'MarkerEdgeColor','r');
    set(mplot(2,options.sides),'MarkerEdgeColor','g');

setappdata(gcf,'selectrode',0);

function ea_view(hobj,ev,commnd)
switch commnd
    case 'p'
        view(0,0);
    case {'x','a'}
        view(180,0);
    case {'y','r'}
        view(90,0);
    case 'l'
        view(270,0);
end

function ea_finish(hobj,ev)
disp('Manual correction done.');


function robotSpace(hobj,ev) % simulates key presses using Java.
import java.awt.Robot;
import java.awt.event.*;
SimKey=Robot;
SimKey.keyPress(KeyEvent.VK_SPACE)
SimKey.keyRelease(KeyEvent.VK_SPACE)

function y = ea_nanmean(varargin)
if nargin==2
    x=varargin{1};
    dim=varargin{2};
elseif nargin==1
x=varargin{1};
    dim=1;
end

N = sum(~isnan(x), dim);
y = ea_nansum(x, dim) ./ N;


function trajectory=ea_prolong_traj(trajectory)
maxv=max([length(trajectory{1}),length(trajectory{2})]);
for side=1:length(trajectory)
    for long=1:maxv-length(trajectory{side})+20
        trajectory{side}(end+1,:)=trajectory{side}(end,:)+(trajectory{side}(end,:)-trajectory{side}(end-1,:));
    end
end


function M = ea_nmi(X,Y)
% function M = MI_GG(X,Y)
% Compute the mutual information of two images: X and Y, having
% integer values.
%
% INPUT:
% X --> first image
% Y --> second image (same size of X)
%
% OUTPUT:
% M --> mutual information of X and Y
%
% Written by GIANGREGORIO Generoso.
% DATE: 04/05/2012
% E-MAIL: ggiangre@unisannio.it
%__________________________________________________________________________

X = double(X);
Y = double(Y);

X_norm = X - min(X(:)) +1;
Y_norm = Y - min(Y(:)) +1;

matAB(:,1) = X_norm(:);
matAB(:,2) = Y_norm(:);
h = accumarray(matAB+1, 1); % joint histogram

hn = h./sum(h(:)); % normalized joint histogram
y_marg=sum(hn,1);
x_marg=sum(hn,2);

Hy = - sum(y_marg.*log2(y_marg + (y_marg == 0))); % Entropy of Y
Hx = - sum(x_marg.*log2(x_marg + (x_marg == 0))); % Entropy of X

arg_xy2 = hn.*(log2(hn+(hn==0)));
h_xy = sum(-arg_xy2(:)); % joint entropy
M = Hx + Hy - h_xy; % mutual information
