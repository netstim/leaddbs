function ea_manualreconstruction(mcfig,patientname,options)
% This is the function that enables the user to manually correct the
% electrode positions that have been reconstructed by the algorithm. A
% handle for the figure, the coordinates of the contacts in millimeter
% notation, optionally manually measured coordinates, the fitted line in
% form of a 1x2 cell each containing a nx3 matrix that describes the line,
% the full path to the coronal nifti file, the name of the patient and the
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
options.loadnativereco = 1; % Load native reco intead of scrf
options.xray=0;
setappdata(mcfig,'options',options);

if ~isfile(options.subj.recon.recon)
    close(mcfig);
    msgbox('Please run pre-Reconstruct module first.');
    return
end

[~,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);

setappdata(mcfig,'origtrajectory',trajectory);
setappdata(mcfig,'manually_corrected',manually_corrected);

% initialize scene
ea_mancor_updatescene([],[],mcfig);

% initialize toolbar
ht=uitoolbar(mcfig);
captions=getappdata(mcfig,'captions');
c_step=2;

prevlead=uipushtool(ht,'CData',ea_get_icn('prevlead'),'TooltipString','Previous Electrode [<]','ClickedCallback',{@sequenceelectrode,mcfig,'prev'});
nextlead=uipushtool(ht,'CData',ea_get_icn('nextlead'),'TooltipString','Next Electrode [>]','ClickedCallback',{@sequenceelectrode,mcfig,'next'});

hf=uitoggletool(ht,'CData',ea_get_icn('hidefid'),'TooltipString','Show/Hide Fiducials','State','off','OnCallback',{@hidefid,mcfig},'OffCallback',{@showfid,mcfig});
setappdata(mcfig,'hf',hf);

minuscontrast=uipushtool(ht,'CData',ea_get_icn('contrastminus'),'TooltipString','Decrease Contrast [C]','ClickedCallback',{@setcontrast,'c',nan,mcfig});
pluscontrast=uipushtool(ht,'CData',ea_get_icn('contrastplus'),'TooltipString','Increase Contrast [V]','ClickedCallback',{@setcontrast,'v',nan,mcfig});
minusbrightness=uipushtool(ht,'CData',ea_get_icn('brightnessminus'),'TooltipString','Decrease Brightness [B]','ClickedCallback',{@setcontrast,'b',nan,mcfig});
plusbrightness=uipushtool(ht,'CData',ea_get_icn('brightnessplus'),'TooltipString','Increase Brightness [N]','ClickedCallback',{@setcontrast,'n',nan,mcfig});

eltog(1)=uitoggletool(ht,'CData',ea_get_icn('el0'),'TooltipString','Select distal contact [0]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
eltog(2)=uitoggletool(ht,'CData',ea_get_icn('el3'),'TooltipString','Select proximal contact [3]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
%eltog(3)=uitoggletool(ht,'CData',ea_get_icn('el4'),'TooltipString','Select Electrode 4 [4]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});
%eltog(4)=uitoggletool(ht,'CData',ea_get_icn('el7'),'TooltipString','Select Electrode 7 [7]','State','off','OnCallback',{@selectelectrode},'OffCallback',{@deselectelectrode});

% if numel(options.sides) == 1    % only one hemisphere
%     switch options.sides
%         case 1  % right hemisphere, disable '4' and '7' buttons
%             set(eltog(3), 'Enable', 'off');
%             set(eltog(4), 'Enable', 'off')
%         case 2  % left hemisphere, disable '0' and '3' buttons
%             set(eltog(1), 'Enable', 'off')
%             set(eltog(2), 'Enable', 'off')
%     end
% end

%xview=uipushtool(ht,'CData',ea_get_icn('elX'),'TooltipString','Set view from X-Direction [X]','ClickedCallback',{@ea_view,'x'});
%yview=uipushtool(ht,'CData',ea_get_icn('elY'),'TooltipString','Set view from Y-Direction [Y]','ClickedCallback',{@ea_view,'y'});

% rotleft=uipushtool(ht,'CData',ea_get_icn('rotleft'),'TooltipString','Rotate Electrode counter-clockwise','ClickedCallback',{@ea_rotate,'cc',mcfig});
% rotright=uipushtool(ht,'CData',ea_get_icn('rotright'),'TooltipString','Rotate Electrode clockwise','ClickedCallback',{@ea_rotate,'c',mcfig});

% rotleftcw=uipushtool(ht,'CData',ea_get_icn('rotleftcw'),'TooltipString','Rotate left lead clockwise','ClickedCallback',{@ea_rotate,'lc',mcfig});
% rotleftccw=uipushtool(ht,'CData',ea_get_icn('rotleftccw'),'TooltipString','Rotate left lead counterclockwise','ClickedCallback',{@ea_rotate,'lcc',mcfig});
autorotation=uipushtool(ht,'CData',ea_get_icn('autorot'),'TooltipString','Detect directional leads','ClickedCallback',{@ea_autorotate,'clockwise',mcfig});
manualrotation=uipushtool(ht,'CData',ea_get_icn('manualrot'),'TooltipString','Detect directional leads','ClickedCallback',{@ea_manualrotate,'clockwise',mcfig});
rotationcw=uipushtool(ht,'CData',ea_get_icn('cw'),'TooltipString','Rotate lead clockwise','ClickedCallback',{@ea_rotate,'clockwise',mcfig});
rotationccw=uipushtool(ht,'CData',ea_get_icn('ccw'),'TooltipString','Rotate lead counterclockwise','ClickedCallback',{@ea_rotate,'counterclockwise',mcfig});

% mni=uitoggletool(ht,'CData',ea_get_icn('mninative'),'TooltipString','Toggle MNI vs. Native space','State','off','OnCallback',{@ea_mancor_updatescene,mcfig,'mni'},'OffCallback',{@ea_mancor_updatescene,mcfig,'native'});

finish_mc=uipushtool(ht,'CData',ea_get_icn('done'),'TooltipString','Finish manual corrections [space]','ClickedCallback',{@robotSpace});

captoggle=uitoggletool(ht,'CData',ea_get_icn('labels'),'TooltipString','Orientation','OnCallback',{@objvisible,captions},'OffCallback',{@objinvisible,captions},'State','on');

xt=uitoggletool(ht,'CData',ea_get_icn('xray'),'TooltipString','Toggle X-Ray Mode','State','off','OnCallback',{@xrayon,mcfig},'OffCallback',{@xrayoff,mcfig});
setappdata(mcfig,'xt',xt);

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
% markers=getappdata(gcf,'markers');
% trajectory=getappdata(gcf,'trajectory');
options=getappdata(mcfig,'origoptions');

options.hybridsave=1;
options.native=1;
options.loadnativereco = 1; % Load native reco intead of scrf
[coords_mm,trajectory,markers]=ea_load_reconstruction(options);

ea_save_reconstruction(coords_mm,trajectory,markers,options.elmodel,1,options);
options=getappdata(mcfig,'origoptions');
try options=rmfield(options,'hybridsave'); end
ea_busyaction('off',mcfig,'reco');
close(mcfig);

if options.autoimprove
    disp('Storing results in template.');
    ea_export_templates(coords_mm{1}(1:4,:),trajectory{1},options.patientname,options,'r')
    ea_export_templates(coords_mm{2}(1:4,:),trajectory{2},options.patientname,options,'l')
    disp('Done.');
end

%% methods dump:
ea_methods(options,...
    ['DBS-Electrodes were manually localized based on post-operative acquisitions using a tool specifically designed for this task (as implemented in Lead-DBS software',...
    '; Horn & Kuehn 2005; SCR_002915; https://www.lead-dbs.org).'],...
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
options.loadnativereco = 1;  % Load native reco intead of scrf
[coords_mm,trajectory,markers,elmodel]=ea_load_reconstruction(options);
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
        %     case 'L0'
        %         markers(2).tail=markers(2).tail+(mLH-markers(2).head);
        %         markers(2).head=mLH;
        %     case 'L3'
        %         markers(2).head=markers(2).head+(mLT-markers(2).tail);
        %         markers(2).tail=mLT;
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
% pause
% commnd=get (gcf, 'CurrentKey');

% get vars
eltog=getappdata(mcfig,'eltog');
elplot=getappdata(mcfig,'elplot');
mplot=getappdata(mcfig,'mplot');
options=getappdata(mcfig,'options');
if ~isfield(options,'visible')
    options.visible=1;
end
if  ~isfile(options.subj.recon.recon)
    close(mcfig);
    return
end

options.loadnativereco = 1;  % Load native reco intead of scrf
[~,~,markers]=ea_load_reconstruction(options);

commnd=event.Character;
switch lower(commnd)
    case ' '
        %markers=getappdata(mcfig,'markers');
        if options.elside==options.sides(end)
            ea_endfcn(mcfig);
            return
        else % move to next electrode in batch
            set(0,'CurrentFigure',mcfig);
            setappdata(mcfig,'selectrode',0);
            [~,ix]=ismember(options.elside,options.sides);
            options.elside=options.sides(ix+1);
            setappdata(mcfig,'options',options);
            ea_mancor_updatescene([],[],mcfig);
        end
    case '.' % center selected electrode - this doesn't seem to work yet.
        selectrode=getappdata(mcfig,'selectrode');
        if selectrode
            
            optoffsets=getappdata(mcfig,'optoffsets');
            [coords_mm,trajectory,markers,elmodel]=ea_load_reconstruction(options);
            movedcoords=moveonecoord(markers,selectrode,optoffsets(selectrode,:),options); % move the correct coord to the correct direction.
            
            set(mplot(1,1),'XData',movedcoords(options.elside).head(1),'YData',movedcoords(options.elside).head(2),'ZData',movedcoords(options.elside).head(3))
            set(mplot(2,1),'XData',movedcoords(options.elside).tail(1),'YData',movedcoords(options.elside).tail(2),'ZData',movedcoords(options.elside).tail(3))
            %            setappdata(mcfig,'markers',markers);
            ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);
            ea_mancor_updatescene([],[],mcfig);
        end
    case '<' % move to previous electrode in batch
        sequenceelectrode([],[],mcfig,'prev');
    case '>' % move to next electrode in batch
        sequenceelectrode([],[],mcfig,'next')
        %     case {'x','a','p','y','l','r'} % view angles.
        %         %markers=getappdata(mcfig,'markers');
        %         ea_view(nan,nan,commnd);
    case 'h'
        showhidefid([],[],mcfig);
    case 'x'
        xrayswitch([],[],mcfig);
    case {'0','3','4','7'}
        switch lower(commnd)
            case {'0','4'}
                selectrode=1;
            case {'3','7'}
                selectrode=2;
        end
        
        if numel(options.sides) == 1	% only one hemisphere
            switch options.sides
                case 1	% right hemisphere, return if '4' or '7' is pressed
                    if selectrode == 3 || selectrode == 4	%
                        return;
                    end
                case 2	% left hemisphere, return if '0' or '3' is pressed
                    if selectrode == 1 || selectrode == 2
                        return;
                    end
            end
        end
        
        oselectrode=getappdata(mcfig,'selectrode');
        if selectrode==oselectrode % toggle had already been clicked -> deselect all.
            % reset all toggletools
            for i=1:2
                set(eltog(i),'State','off');
            end
            setappdata(mcfig,'selectrode',0);
            ea_mancor_updatescene([],[],mcfig);
        else
            % clear all toggletools.
            for i=1:2
                set(eltog(i),'State','off');
            end
            
            % set the correct toggletool again.
            set(eltog(selectrode),'State','on');
            
            % store selected electrode in appdata.
            setappdata(mcfig,'selectrode',selectrode);
            ea_mancor_updatescene([],[],mcfig);
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
        ea_mancor_updatescene([],[],mcfig);
        
        % reload the parameters
        [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
    otherwise % arrow keys, plus, minus
        if ismember(event.Key,{'rightarrow','leftarrow','uparrow','downarrow'}) || ismember(event.Character,{'+','-','*','_'})
            selectrode=getappdata(mcfig,'selectrode');
            if ~selectrode % no electrode is highlighted, move electrodes alongside trajectory or increase/decrease spacing.
                [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
                
                markers=ea_correctcoords(markers,trajectory,event,options);
                if isfield(options,'hybridsave')
                    options=rmfield(options,'hybridsave');
                end
                ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);
                
                ea_mancor_updatescene([],[],mcfig);
                [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
                
            else % electrode is highlighted. Move in xy dirs.
                [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
                movedcoords=moveonecoord(markers,selectrode,event,options); % move the correct coord to the correct direction.
                
                set(mplot(1,1),'XData',movedcoords(options.elside).head(1),'YData',movedcoords(options.elside).head(2),'ZData',movedcoords(options.elside).head(3))
                set(mplot(2,1),'XData',movedcoords(options.elside).tail(1),'YData',movedcoords(options.elside).tail(2),'ZData',movedcoords(options.elside).tail(3))
                
                ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,1,options);
                
                ea_mancor_updatescene([],[],mcfig);
                
                [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
            end
        end
end

cnt=1;
options=getappdata(mcfig,'options');
coords_mm=ea_resolvecoords(markers,options);

set(mplot(1,1),'XData',markers(options.elside).head(1),'YData',markers(options.elside).head(2),'ZData',markers(options.elside).head(3));
set(mplot(2,1),'XData',markers(options.elside).tail(1),'YData',markers(options.elside).tail(2),'ZData',markers(options.elside).tail(3));
for el=1:length(elplot)/2
    set(elplot(cnt),'XData',coords_mm{options.elside}(el,1),'YData',coords_mm{options.elside}(el,2),'ZData',coords_mm{options.elside}(el,3))
    cnt=cnt+1;
end
refreshdata(elplot,'caller')
drawnow


function showhidefid(m,n,mcfig)
options=getappdata(mcfig,'options');
if options.visible
    hidefid([],[],mcfig);
else
    showfid([],[],mcfig);
end


function hidefid(m,n,mcfig)
options=getappdata(mcfig,'options');
options.visible=0;
hf=getappdata(mcfig,'hf');
setappdata(mcfig,'options',options);
ea_mancor_updatescene([],[],mcfig);
set(hf,'State','on');


function showfid(m,n,mcfig)
options=getappdata(mcfig,'options');
options.visible=1;
hf=getappdata(mcfig,'hf');
setappdata(mcfig,'options',options);
ea_mancor_updatescene([],[],mcfig);
set(hf,'State','off');


function xrayswitch(hf,t,mcfig)
options=getappdata(mcfig,'options');
if options.xray
    xrayoff(hf,t,mcfig);
else
    xrayon(hf,t,mcfig);
end


function xrayoff(hf,t,mcfig)
options=getappdata(mcfig,'options');
options.xray=0;
xt=getappdata(mcfig,'xt');
setappdata(mcfig,'options',options);
ea_mancor_updatescene([],[],mcfig);
set(xt,'State','off');


function xrayon(hf,t,mcfig)
options=getappdata(mcfig,'options');
options.xray=1;
xt=getappdata(mcfig,'xt');
setappdata(mcfig,'options',options);
ea_mancor_updatescene([],[],mcfig);
set(xt,'State','on');


function sequenceelectrode(m,n,mcfig,what)
options=getappdata(mcfig,'options');
switch what
    case 'next'
        [~,ix]=ismember(options.elside,options.sides);
        if ix<length(options.sides)
            setappdata(mcfig,'selectrode',0);
            options.elside=options.sides(ix+1);
            setappdata(mcfig,'options',options);
            ea_mancor_updatescene([],[],mcfig);
        end
    case 'prev'
        [~,ix]=ismember(options.elside,options.sides);
        if ix>1
            setappdata(mcfig,'selectrode',0);
            options.elside=options.sides(ix-1);
            setappdata(mcfig,'options',options);
            ea_mancor_updatescene([],[],mcfig);
        end
end


function sp=getsuplots(sides)
if isequal(sides,[1:2])
    sp=1:4;
elseif isequal(sides,1)
    sp=1:2;
elseif isequal(sides,2)
    sp=3:4;
end


function markers=moveonecoord(markers,selectrode,command,options)

grone=[0.02,0.1,1]; % step-sizes

movedcoords=[markers(options.elside).head;markers(options.elside).tail];

if isnumeric(command)
    movedcoords(selectrode,:)=movedcoords(selectrode,:)+[command,0];
else
    switch command.Key
        case 'leftarrow'
            movedcoords(selectrode,1)=movedcoords(selectrode,1)-grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier));
        case 'rightarrow'
            movedcoords(selectrode,1)=movedcoords(selectrode,1)+grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier));
        case 'uparrow'
            movedcoords(selectrode,2)=movedcoords(selectrode,2)+grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier));
        case 'downarrow'
            movedcoords(selectrode,2)=movedcoords(selectrode,2)-grone(2+ismember('shift',command.Modifier)-ismember('alt',command.Modifier));
    end
end
markers(options.elside).head=movedcoords(1,:);
markers(options.elside).tail=movedcoords(2,:);


function objvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');
%disp([atls,'visible clicked']);


function objinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');
%disp([atls,'invisible clicked']);


function ea_rotate(hobj,ev,ccw,mcfig)
rotation=getappdata(gcf,'rotation'); % rotation angle in degrees
options = getappdata(gcf,'options');
if isempty(rotation)
    rotation{1} = 0;
    rotation{2} = 0;
end
sideact = options.elside;
switch ccw
    case 'clockwise'
        rotation{sideact} = rotation{sideact} + 15;
    case 'counterclockwise'
        rotation{sideact} = rotation{sideact} - 15;
end
setappdata(gcf,'rotation',rotation);
ea_mancor_updatescene([],[],mcfig);


function ea_autorotate(hobj,ev,ccw,mcfig)
options = getappdata(gcf,'options');
rotation=getappdata(gcf,'rotation'); % rotation angle in degrees
orientation = ea_diode_main(options);
if ~isempty(orientation)
    rotation{options.elside} = orientation;
end
figure(mcfig);
setappdata(gcf,'rotation',rotation);
ea_mancor_updatescene([],[],mcfig);

function ea_manualrotate(hobj,ev,ccw,mcfig)
options = getappdata(gcf,'options');
rotation=getappdata(gcf,'rotation'); % rotation angle in degrees
orientation = ea_diode_manual_main(options);
if ~isempty(orientation)
    rotation{options.elside} = orientation;
end
figure(mcfig);
setappdata(gcf,'rotation',rotation);
ea_mancor_updatescene([],[],mcfig);

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
if contrast<0.1
    contrast=0.1;
end
setappdata(mcfig,'contrast',contrast);
setappdata(mcfig,'offset',offset);
%disp(['Contrast: ',num2str(contrast),', Offset: ',num2str(offset),'.']);
ea_mancor_updatescene([],[],mcfig);


function selectelectrode(hobj,ev)
% reset all toggletools
warning ('off','all');
eltog=getappdata(gcf,'eltog');
elplot=getappdata(gcf,'elplot');
mplot=getappdata(gcf,'mplot');
options=getappdata(gcf,'options');
set(mplot(1,1),'MarkerEdgeColor','r');
set(mplot(2,1),'MarkerEdgeColor','g');

for i=getsuplots(1)
    set(eltog(i),'State','off');
    if eltog(i)==hobj
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
set(mplot(1,1),'MarkerEdgeColor','r');
set(mplot(2,1),'MarkerEdgeColor','g');

setappdata(gcf,'selectrode',0);


function robotSpace(hobj,ev) % simulates key presses using Java.
import java.awt.Robot;
import java.awt.event.*;
SimKey=Robot;
SimKey.keyPress(KeyEvent.VK_SPACE)
SimKey.keyRelease(KeyEvent.VK_SPACE)
