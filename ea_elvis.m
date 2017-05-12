function resultfig=ea_elvis(varargin)
% This function is the main 3D-visualization function of lead DBS. It
% reads all atlases found in the lead_root/atlases folder, calculates a
% convex hull around the nonzero area and renders it as 3D surfaces.
% Resulting figures are interactive and objects can be hidden. This
% functionality is preserved as long as the files are saved as .fig Matlab
% figures and the functions remain on the path when figures are reopened
% lateron. Basically, this function only needs a valid options struct as input.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
%
% Andreas Horn


% Initialize inputs
options=varargin{1};
if nargin>2
    stimparams=varargin{3};
else
    stimparams=nan;
end
if nargin==4
    fiberthresh=varargin{4};
else

    fiberthresh=options.fiberthresh;
end

% Initialize figure

resultfig=figure('name',[options.patientname,': Electrode-Scene'],'color','k','numbertitle','off','CloseRequestFcn',@closesattelites,'visible',options.d3.verbose,'KeyPressFcn',@ea_keypress,'KeyReleaseFcn',@ea_keyrelease);
setappdata(resultfig,'options',options);
set(resultfig,'toolbar','none');
ssz=get(0,'Screensize');
ssz(1:2)=ssz(1:2)+50;
ssz(3:4)=ssz(3:4)-200;
set(resultfig, 'Position', ssz); % Maximize figure.

% initialize some ui elements
ht=uitoolbar(resultfig);

% add custom rotator:
uibjs.rotate3dtog=uitoggletool(ht,'CData',ea_get_icn('rotate',options),'TooltipString','Rotate 3D','OnCallback',{@ea_rotate,'on'},'OffCallback',{@ea_rotate,'off'},'State','off');
uibjs.magnifyplus=uitoggletool(ht,'CData',ea_get_icn('magnplus',options),'TooltipString','Zoom In','OnCallback',{@ea_zoomin,'on'},'OffCallback',{@ea_zoomin,'off'},'State','off');
uibjs.magnifyminus=uitoggletool(ht,'CData',ea_get_icn('magnminus',options),'TooltipString','Zoom Out','OnCallback',{@ea_zoomout,'on'},'OffCallback',{@ea_zoomout,'off'},'State','off');
uibjs.handtog=uitoggletool(ht,'CData',ea_get_icn('hand',options),'TooltipString','Pan Scene','OnCallback',{@ea_pan,'on'},'OffCallback',{@ea_pan,'off'},'State','off');
setappdata(resultfig,'uibjs',uibjs);


mh = uimenu(resultfig,'Label','Add Objects');
fh1 = uimenu(mh,'Label','Open Tract','Callback',{@ea_addobj,resultfig,'tract',options});
fh2 = uimenu(mh,'Label','Open ROI','Callback',{@ea_addobj,resultfig,'roi',options});
fh3 = uimenu(mh,'Label','Show tracts weighted by activation map','Callback',{@ea_addobj,resultfig,'tractmap',options});

% Set some visualization parameters
set(resultfig,'Renderer','opengl')
axis off
%set(gca,'DrawMode','fast')
set(resultfig, 'InvertHardCopy', 'off');
%set(resultfig,'visible','off');
set(resultfig,'Clipping','off');
set(gca,'NextPlot','replacechildren');
%set(gca,'Erasemode','none');

% Get paramters (Coordinates and fitted line).
figtitle=get(gcf,'Name');
set(gcf,'Name',[figtitle,'...building...']);
axis equal
axis fill

colormap('gray')



%% Patient specific part (skipped if no patient is selected or no reco available):
if ~strcmp(options.patientname,'No Patient Selected') % if not initialize empty viewer
    if exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file') || nargin>1;
        if nargin>1
            multiplemode=1;
            elstruct=varargin{2};
            
            if options.d3.mirrorsides
               elstruct=ea_mirrorsides(elstruct); 
               options.d3.isomatrix=ea_mirrorsides(options.d3.isomatrix);
            end
            
        else
            multiplemode=0;
            options.loadrecoforviz=1; % add flag to load scrf entry if in native mode.
            [coords_mm,trajectory,markers]=ea_load_reconstruction(options);

            elstruct(1).coords_mm=coords_mm;
            elstruct(1).coords_mm=ea_resolvecoords(markers,options);
            elstruct(1).trajectory=trajectory;
            elstruct(1).name=options.patientname;
            elstruct(1).markers=markers;
            clear coords_mm trajectory
        end

        % show electrodes..
        for pt=1:length(elstruct)
            try
                [el_render(pt).el_render,el_label(:,pt)]=ea_showelectrode(resultfig,elstruct(pt),pt,options);
            catch
                ea_error(['Couldn''t visualize electrode from patient ',num2str(pt),'.']);
            end
            if options.d3.elrendering==1 % export vizstruct for lateron export to JSON file / Brainbrowser.
               % this part for brainbrowser support.
               vizstruct=struct('faces',[],'vertices',[],'colors',[]);

               cnt=1;
                for side=1:length(options.sides)
                    extract=1:length(el_render(pt).el_render{side});
                    for ex=extract
                        tp=el_render(pt).el_render{side}(ex);

                        try % works only in ML 2015
                            tr=triangulation(get(el_render(pt).el_render{side}(ex),'Faces'),get(el_render(pt).el_render{side}(ex),'Vertices'));
                            vizstruct(cnt).normals = vertexNormal(tr);
                        catch % workaround for older versions..
                            vizstruct(cnt).normals=get(tp,'VertexNormals')*-1;
                        end

                        vizstruct(cnt).faces=get(tp,'Faces');
                        vizstruct(cnt).vertices=get(tp,'Vertices');
                        scolor=get(el_render(pt).el_render{side}(ex),'FaceVertexCData');
                        vizstruct(cnt).colors=scolor;
                        %vizstruct(cnt).colors=repmat([squeeze(scolor(1,1,:))',0.7],length(vizstruct(cnt).faces),1);
                        vizstruct(cnt).name='';
                        cnt=cnt+1;
                    end
                end
            end
        end

        setappdata(resultfig,'el_render',el_render);
        % add handles to buttons. Can't be combined with the above loop since all
        % handles need to be set for the buttons to work properly (if alt is
        % pressed, all electrodes are made visible/invisible).
        drawnow

        try
            set(el_label,'Visible','off');
            ellabeltog=uitoggletool(ht,'CData',ea_get_icn('labels',options),'TooltipString','Electrode labels','OnCallback',{@objvisible,el_label},'OffCallback',{@objinvisible,el_label},'State','off');
        end

        cnt=1;
        for pt=1:length(elstruct)
                try
                    if multiplemode
                        caption{1}=[elstruct(pt).name,'_Left'];         caption{2}=[elstruct(pt).name,'_Right'];
                    else
                        caption{1}='Electrode_Left'; caption{2}='Electrode_Right';
                    end
                    eltog(cnt)=uitoggletool(ht,'CData',ea_get_icn('electrode',options),'TooltipString',caption{1},'OnCallback',{@elvisible,el_render,pt,2,'on',options},'OffCallback',{@elvisible,el_render,pt,2,'off',options},'State','on');
                    eltog(cnt+1)=uitoggletool(ht,'CData',ea_get_icn('electrode',options),'TooltipString',caption{2},'OnCallback',{@elvisible,el_render,pt,1,'on',options},'OffCallback',{@elvisible,el_render,pt,1,'off',options},'State','on');
                    
                    if exist([options.uipatdirs{pt} '/cortex/CortElecs.mat'],'file')
                    vars = whos('-file',[options.uipatdirs{pt} '/cortex/CortElecs.mat']);
                    CortElecs = load([options.uipatdirs{pt} '/cortex/CortElecs.mat']);    
                        if ismember('Left',{vars.name})
                            hold on; plot3(CortElecs.Left(:,1),CortElecs.Left(:,2),CortElecs.Left(:,3),'.','color','r','markersize',10)
                            ctxeltog(cnt)=uitoggletool(ht,'CData',ea_get_icn('cortical_strip',options),'TooltipString',['Cortical_' caption{1}],'OnCallback',{@ctxelvisible,el_renderstrip,pt,2,'on',options},'OffCallback',{@elvisible,el_render,pt,2,'off',options},'State','on');
                        end
                        if ismember('Right',{vars.name})
                            ctxeltog(cnt)=uitoggletool(ht,'CData',ea_get_icn('cortical_strip',options),'TooltipString',['Cortical_' caption{1}],'OnCallback',{@ctxelvisible,el_renderstrip,pt,2,'on',options},'OffCallback',{@elvisible,el_render,pt,2,'off',options},'State','on');
                        end
                    end
                    cnt=cnt+2;
                end
        end
        setappdata(resultfig,'eltog',eltog);

        clear cnt

        % Initialize Stimulation-Button

        stimbutton=uipushtool(ht,'CData',ea_get_icn('stimulation',options),'TooltipString','Stimulation Control Figure','ClickedCallback',{@openstimviewer,elstruct,resultfig,options});

    else
        options.writeoutstats=0; % if no electrodes are there, stats can't be written.
        elstruct=struct;
    end
    
else
    options.writeoutstats=0; % if no electrodes are there, stats can't be written.
    elstruct=struct;
end



% Initialize Sliceview-Button

slicebutton=uipushtool(ht,'CData',ea_get_icn('slices',options),'TooltipString','Slice Control Figure','ClickedCallback',{@opensliceviewer,resultfig,options});

if options.prefs.env.dev;
% Initialize MER-Button
merbutton=uipushtool(ht,'CData',ea_get_icn('mer',options),'TooltipString','MER Control Figure','ClickedCallback',{@openmerviewer,resultfig,options});
end
% Initialize Convis-Button
convisbutton=uipushtool(ht,'CData',ea_get_icn('connectome',options),'TooltipString','Connectivity Visualization','ClickedCallback',{@openconnectomeviewer,resultfig,options});

% Initialize FS Cortex-Button
corticalbutton=uipushtool(ht,'CData',ea_get_icn('cortex',options),'TooltipString','Cortical Reconstruction Visualization','ClickedCallback',{@opencortexviewer,resultfig,options});

% Initialize Cortical Strip-Button
% cortelsbutton=uipushtool(ht,'CData',ea_get_icn('cortical_strip',options),'TooltipString','Cortical Reconstruction Visualization','ClickedCallback',{@opencortelsviewer,resultfig,options});

% Show atlas data
if options.d3.writeatlases
    [atlases,colorbuttons,atlassurfs]=ea_showatlas(resultfig,elstruct,options);
    
    %if length(atlases.names)>6 % only open up for long atlas lists by default.
    ea_openatlascontrol([],[],atlases,resultfig,options);
    %end

    if options.d3.elrendering==1 % export vizstruct for lateron export to JSON file / Brainbrowser.
        try % see if electrode has been defined.
        cnt=length(vizstruct);
        catch
            cnt=0;
        end
        % export vizstruct
        try
            for side=1:length(options.sides)
                for atl=1:length(atlases.fv)

                    if isfield(atlases.fv{atl,side},'faces')
                        vizstruct(cnt+1).faces=atlases.fv{atl,side}.faces;
                        vizstruct(cnt+1).vertices=atlases.fv{atl,side}.vertices;
                        vizstruct(cnt+1).normals=atlases.normals{atl,side};
                        vizstruct(cnt+1).colors=[squeeze(ind2rgb(round(atlases.cdat{atl,side}),atlases.colormap)),repmat(0.7,size(atlases.normals{atl,side},1),1)];
                        cnt=cnt+1;
                    end
                end
            end
        end
    end
end

% Show isomatrix data

if options.d3.showisovolume
    allisomatrices=options.d3.isomatrix;
    allisonames=options.d3.isomatrix_name;
    for reg=1:length(allisomatrices)
        options.d3.isomatrix=allisomatrices{reg};
        options.d3.isomatrix_name=allisonames{reg};
        

        ea_showisovolume(resultfig,elstruct,options);
    end
end

if isfield(options.d3,'expdf');
    if options.d3.expdf
        %cd([options.root,options.patientname]);
        fig2pdf3d(gca,[options.root,options.patientname,filesep,'Lead-DBS_Electrode_Localization'],options);
        close(resultfig);
        return
    end
end
%% End of patient-specific part.

% Initialize a draggable lightbulb
hold on
ea_show_light(resultfig,1);
% set(lightbulb, 'Visible', 'off');

lightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('lightbulb',options),'TooltipString','Lightbulb','OnCallback',{@objvisible,getappdata(resultfig,'cam_lamp')},'OffCallback',{@objinvisible,getappdata(resultfig,'cam_lamp')},'State','on');
clightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('clightbulb',options),'TooltipString','Lightbulb','OnCallback',{@objvisible,getappdata(resultfig,'ceiling_lamp')},'OffCallback',{@objinvisible,getappdata(resultfig,'ceiling_lamp')},'State','on');
llightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('llightbulb',options),'TooltipString','Lightbulb','OnCallback',{@objvisible,getappdata(resultfig,'left_lamp')},'OffCallback',{@objinvisible,getappdata(resultfig,'left_lamp')},'State','on');
rlightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('rlightbulb',options),'TooltipString','Lightbulb','OnCallback',{@objvisible,getappdata(resultfig,'right_lamp')},'OffCallback',{@objinvisible,getappdata(resultfig,'right_lamp')},'State','on');


% Initialize HD-Export button

hdsavebutton=uipushtool(ht,'CData',ea_get_icn('save',options),'TooltipString','Save Scene','ClickedCallback',@export_hd);
dofsavebutton=uipushtool(ht,'CData',ea_get_icn('save_depth',options),'TooltipString','Save Scene with depth of field','ClickedCallback',{@ea_export_depth_of_field,resultfig});


% Initialize Video-Export button

videoexportbutton=uipushtool(ht,'CData',ea_get_icn('video',options),'TooltipString','Save video','ClickedCallback',{@export_video,options});


% Init hard_electrode_view button
if options.modality==2
electrodesegmentbutton=uitoggletool(ht,'CData',ea_get_icn('electrode_segment',options),'TooltipString','Auto-Segment electrode from postoperative acquisition','OnCallback',{@ea_segment_electrode,options,resultfig,'on'},'OffCallback',{@ea_segment_electrode,options,resultfig,'off'},'State','off');
end
% Initialize Export to Lead-Server button

lsbutton=uipushtool(ht,'CData',ea_get_icn('server',options),'TooltipString','Export to Server','ClickedCallback',{@ea_export_server,options});

hold off

set(0,'CurrentFigure',resultfig);

set(resultfig,'Renderer','OpenGL')
axis off
set(resultfig,'color','k');

axis vis3d
axis equal
set(resultfig,'Name',figtitle);
set(0,'CurrentFigure',resultfig);
ax=resultfig.CurrentAxes;
set(ax,'XLimMode','auto');
set(ax,'YLimMode','auto');
set(ax,'ZLimMode','auto');

% set(ax,'XLim',[-140 140]);
% set(ax,'YLim',[-140 140]);
% set(ax,'ZLim',[-140 140]);
% zoom(3)
% ax=gca;
% set(ax,'XLim',[-140 140]);
% set(ax,'YLim',[-140 140]);
% set(ax,'ZLim',[-140 140]);
view(142,13.6)
zoom(1.5)
%set(resultfig,'visible','on');
opensliceviewer([],[],resultfig,options);
if options.d3.elrendering==1 % export vizstruct for lateron export to JSON file / Brainbrowser.
    try
        % store json in figure file
        bbstruct=ea_viz2brainbrowser(vizstruct);
        setappdata(resultfig,'bbstruct',bbstruct);
    end
    if options.prefs.ls.autosave
        ea_export_server([],[],options);
    end
end
setappdata(resultfig,'elstruct',elstruct);


function opensliceviewer(hobj,ev,resultfig,options)
awin=ea_anatomycontrol(resultfig,options);
setappdata(resultfig,'awin',awin);
try WinOnTop(awin,true); end

function openconnectomeviewer(hobj,ev,resultfig,options)
conwin=ea_convis(gcf,options);
setappdata(resultfig,'conwin',conwin);
try WinOnTop(conwin,true); end


function openstimviewer(hobj,ev,elstruct,resultfig,options)
stimwin=ea_stimparams(elstruct,gcf,options);
setappdata(resultfig,'stimwin',stimwin);
try WinOnTop(stimwin,true); end

function openmerviewer(hobj,ev,resultfig,options)
merwin=ea_mercontrol(resultfig,options);
setappdata(resultfig,'merwin',merwin);
try WinOnTop(merwin,true); end

function opencortexviewer(hobj,ev,resultfig,options)
cortex=ea_showcortex(resultfig,options);
setappdata(resultfig,'cortex',cortex);
% reload slice viewer to update opacity control
awin=ea_anatomycontrol(resultfig,options);
setappdata(resultfig,'awin',awin);
try WinOnTop(awin,true); end



function closesattelites(src,evnt)
stimwin=getappdata(gcf,'stimwin');
try
    close(stimwin)
end
awin=getappdata(gcf,'awin');
try
    close(awin)
end
aswin=getappdata(gcf,'aswin');
try
    close(aswin)
end
conwin=getappdata(gcf,'conwin');
try
    close(conwin)
end
merwin=getappdata(gcf,'merwin');
try
    close(merwin)
end
delete(gcf)


function export_video(hobj,ev,options)

%% Set up recording parameters (optional), and record
[FileName,PathName] = uiputfile('LEAD_Scene.mp4','Save file name for video');
ea_CaptureFigVid(options.prefs.video.path, [PathName,FileName],options.prefs.video.opts);



function export_hd(hobj,ev)

[FileName,PathName] = uiputfile('LEAD_Scene.png','Save file name');
if FileName
set(gcf, 'Color', [1,1,1]);
[~, cdata] = ea_myaa([4, 2]);

imwrite(cdata, [PathName,FileName], 'png');
end


function objvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');



function objinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');


function elvisible(hobj,ev,atls,pt,side,onoff,options)

if(getappdata(gcf,'altpressed'))

    eltog=getappdata(gcf,'eltog');
    set(eltog,'State',onoff);
    for el=1:length(atls)
        for side=1:length(options.sides)
           try
               set(atls(el).el_render{side}, 'Visible', onoff);
           end
        end
    end
else
set(atls(pt).el_render{side}, 'Visible', onoff);
end

function ctxelvisible(hobj,ev,atls,pt,side,onoff,options)

if(getappdata(gcf,'altpressed'))

    eltog=getappdata(gcf,'eltog');
    set(eltog,'State',onoff);
    for el=1:length(atls)
        for side=1:length(options.sides)
           try
               set(atls(el).el_render{side}, 'Visible', onoff);
           end
        end
    end
else
set(atls(pt).el_render{side}, 'Visible', onoff);
end


function res=ea_if(condition)
res=0;
if ~isempty(condition)
    if condition
        res=1;
    end
end


function ea_keypress(resultfig,event)
% this is the main keypress function for the resultfigure. Add event
% listeners here.
if ismember('alt',event.Modifier)
     setappdata(resultfig,'altpressed',1);
%    disp('Altpressed');
elseif ismember('shift',event.Modifier)
     setappdata(resultfig,'shiftpressed',1);
end

try
    merwin = getappdata(resultfig,'merwin');
    options = getappdata(merwin,'options');
    %merstruct = getappdata(resultfig,'merstruct');
    merhandles = getappdata(resultfig,'merhandles');
    mermarkers = getappdata(resultfig,'mermarkers');
    keymer = getappdata(resultfig,'keymer');
catch
    merwin=[];
end
if ~isempty(merwin) && isvalid(merwin)
    commnd=event.Key; % event.Character;
    sphere = load([options.earoot,'helpers',filesep,'sphere.mat']);
    n = length(mermarkers);
    sSize = 0.5;
    % CData = parula; CData = repmat(CData(1:length(sphere.x),:),[1 size(sphere.x,2)/3]);
    % colormap = repmat([1 1 0],[size(sphere.x,1) size(sphere.x,2)/3]);
    % tmp = parula;
    % colormap = repmat(tmp(end:-2:1,:),[4 1]); clear tmp
    %colormap = [1 1 0; 0 0 1];
    
    if isempty(keymer)
        return
    elseif sum(double(~cellfun(@isempty,strfind({'space','m','l','t','b','s','n'},event.Key))))>0
        trajectory = [get(merhandles.(keymer(4:end)),'XData')',get(merhandles.(keymer(4:end)),'YData')',get(merhandles.(keymer(4:end)),'ZData')'];
        if n>0 && isequal(trajectory(1,:),mermarkers(n).coords_mm) && ~strcmp(event.Key,'s') && ~strcmp(event.Key,'n')
            fprintf('Location along %s %s tract already marked: [%f,%f,%f].\n',keymer(strfind(keymer,'_')+1:end),keymer(4:strfind(keymer,'_')-1),trajectory(1,:))
            return
        end
    end   
    
    if sum(double(~cellfun(@isempty,strfind({'space','m','l','t','b','s','n'},event.Key))))>0
        sphere.x = sphere.x*sSize+trajectory(1,1);
        sphere.y = sphere.y*sSize+trajectory(1,2);
        sphere.z = sphere.z*sSize+trajectory(1,3);
        mermarkers(n+1).side = keymer(regexp(keymer,'_')+1:end);
        mermarkers(n+1).tract = keymer(4:regexp(keymer,'_')-1);
        mermarkers(n+1).depth = str2double(getfield(getfield(getappdata(merwin,'UsedByGUIData_m'),['pos' keymer(4:end)]),'String'));
        mermarkers(n+1).coords_mm = trajectory(1,:);
        mermarkers(n+1).dat.implantedtract = getfield(getfield(getappdata(merwin,'UsedByGUIData_m'),['popupimplantedtract' keymer(regexp(keymer,'_'):end)]),'String');
        mermarkers(n+1).dat.implantedtract = mermarkers(n+1).dat.implantedtract{getfield(getfield(getappdata(merwin,'UsedByGUIData_m'),['popupimplantedtract' keymer(regexp(keymer,'_'):end)]),'Value')};
        mermarkers(n+1).dat.leaddepth = str2double(getfield(getfield(getappdata(merwin,'UsedByGUIData_m'),['editimplanteddepth' keymer(regexp(keymer,'_'):end)]),'String'));
        mermarkers(n+1).dat.offset = event.Key;
        mermarkers(n+1).dat.key = event.Key;
        mermarkers(n+1).tag.depth = getfield(getfield(getappdata(merwin,'UsedByGUIData_m'),['pos' keymer(4:end)]),'String');
        mermarkers(n+1).tag.visible = options.prefs.mer.tag.visible;
        [mermarkers(n+1).tag.handle,mermarkers(n+1).tag.string] = ea_setmertag(mermarkers(n+1).tag,keymer,trajectory(1,:));

        % Reserved keys: {'space','m','l','t','b','s','n'}
        % 'space' = Generic; 'm' = MER; 'l' = LFP; 't' = Top; 'b' = Bottom
        switch lower(commnd)
            case 'space'
                mermarkers(n+1).notes;
                mermarkers(n+1).handle = surf(sphere.x,sphere.y,sphere.z,...
                    'FaceColor',[0.5 0.5 0],'EdgeColor','none',...
                    'FaceAlpha',0.7,'tag','Generic');
                mermarkers(n+1).markertype = 'Generic';
            case 'm'
                mermarkers(n+1).handle = surf(sphere.x,sphere.y,sphere.z,...
                    'FaceColor',[0.5 0 0],'EdgeColor','none',...
                    'FaceAlpha',0.7,'tag','MER');
                mermarkers(n+1).markertype = 'MER recording';
                mermarkers(n+1).session = char(inputdlg('Enter Session'));
            case 'l'
                mermarkers(n+1).handle = surf(sphere.x,sphere.y,sphere.z,...
                    'FaceColor',[0 0.5 0],'EdgeColor','none',...
                    'FaceAlpha',0.7,'tag','LFP');
                mermarkers(n+1).markertype = 'LFP recording';
                mermarkers(n+1).session = char(inputdlg('Enter Session'));
            case 't'
                mermarkers(n+1).handle = surf(sphere.x,sphere.y,sphere.z,...
                    'FaceColor',[0 0 0.5],'EdgeColor','none',...
                    'FaceAlpha',0.7,'tag','Top');
                mermarkers(n+1).markertype = 'Top border';
            case 'b'
                mermarkers(n+1).handle = surf(sphere.x,sphere.y,sphere.z,...
                    'FaceColor',[0 0 0.5],'EdgeColor','none',...
                    'FaceAlpha',0.7,'tag','Bottom');
                mermarkers(n+1).markertype = 'Bottom border';
            case 's'
                mermarkers(n).session = char(inputdlg('Enter Session'));
            case 'n'
                mermarkers(n).notes = char(inputdlg('Enter Notes'));
        end

        setappdata(resultfig,'mermarkers',mermarkers);
        ea_updatemercontrol(keymer,getappdata(merwin,'UsedByGUIData_m'),mermarkers,resultfig,options)
   
    else
        switch commnd
            case {'uparrow','leftarrow'}
            if isempty(keymer); return
            else
              trajectory = [merhandles.(keymer(4:end)).XData',merhandles.(keymer(4:end)).YData',merhandles.(keymer(4:end)).ZData'];            
            end
            if isempty(event.Modifier) || ~ismember(event.Modifier,{'shift','alt'})
                d=0.25; % step size
            elseif ismember(event.Modifier,'shift')
                d=0.75;    % large step
            elseif ismember(event.Modifier,'alt')
                d=0.05;    % small step size
            end
            newtrajectory = ea_getmertrajectory(trajectory,d,options.prefs.mer.length,50);
            ea_updatemertrajectory(getappdata(merwin,'UsedByGUIData_m'),newtrajectory,d,keymer(4:end))
        case {'downarrow','rightarrow'}
            if isempty(keymer); return
            else
              trajectory = [merhandles.(keymer(4:end)).XData',merhandles.(keymer(4:end)).YData',merhandles.(keymer(4:end)).ZData'];
            end
            if isempty(event.Modifier) || ~ismember(event.Modifier,{'shift','alt'})
                d=-0.25; % movement distance
            elseif ismember(event.Modifier,'shift')
                d=-0.75; % large step
            elseif ismember(event.Modifier,'alt')
                d=-0.05;
            end
            newtrajectory = ea_getmertrajectory(trajectory,d,options.prefs.mer.length,50);
            ea_updatemertrajectory(getappdata(merwin,'UsedByGUIData_m'),newtrajectory,d,keymer(4:end))
        end
    end
end
% commnd=event.Character;
% switch lower(commnd)
%     case ' '
%     case {'x','a','p','y','l','r'} % view angles.
%     case {'0','3','4','7'}
%     case {'c','v','b','n'}
%     otherwise % arrow keys, plus, minus
% end


function ea_keyrelease(resultfig,event)
setappdata(resultfig,'altpressed',0);
%disp('Altunpressed');


function [varargout] = ea_myaa(varargin)
% This function has been slightly modified for export use in eAuto-DBS.
% Copyright (c) 2009, Anders Brun
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

%   See also PUBLISH, PRINT
%
%   Version 1.1, 2008-08-21
%   Version 1.0, 2008-08-05
%
%   Author: Anders Brun
%           anders@cb.uu.se
%

%% Force drawing of graphics
drawnow;

%% Find out about the current DPI...
screen_DPI = get(0,'ScreenPixelsPerInch');

%% Determine the best choice of convolver.
% If IPPL is available, imfilter is much faster. Otherwise it does not
% matter too much.
try
    if ippl()
        myconv = @imfilter;
    else
        myconv = @conv2;
    end
catch
    myconv = @conv2;
end

%% Set default options and interpret arguments
if isempty(varargin)
    self.K = [4 4];
    try
        imfilter(zeros(2,2),zeros(2,2));
        self.aamethod = 'imresize';
    catch
        self.aamethod = 'standard';
    end
    self.figmode = 'figure';
elseif strcmp(varargin{1},'publish')
    self.K = [4 4];
    self.aamethod = 'noshrink';
    self.figmode = 'publish';
elseif strcmp(varargin{1},'update')
    self = get(gcf,'UserData');
    figure(self.source_fig);
    drawnow;
    self.figmode = 'update';
elseif strcmp(varargin{1},'lazyupdate')
    self = get(gcf,'UserData');
    self.figmode = 'lazyupdate';
elseif length(varargin) == 1
    self.K = varargin{1};
    if length(self.K) == 1
        self.K = [self.K self.K];
    end
    if self.K(1) > 16
        ea_error('To avoid excessive use of memory, K has been limited to max 16. Change the code to fix this on your own risk.');
    end
    try
        imfilter(zeros(2,2),zeros(2,2));
        self.aamethod = 'imresize';
    catch
        self.aamethod = 'standard';
    end
    self.figmode = 'figure';
elseif length(varargin) == 2
    self.K = varargin{1};
    self.aamethod = varargin{2};
    self.figmode = 'figure';
elseif length(varargin) == 3
    self.K = varargin{1};
    self.aamethod = varargin{2};
    self.figmode = varargin{3};
    if strcmp(self.figmode,'publish') && ~strcmp(varargin{2},'noshrink')
        printf('\nThe AAMETHOD was not set to ''noshrink'': Fixed.\n\n');
        self.aamethod = 'noshrink';
    end
else
    ea_error('Wrong syntax, run: help myaa');
end

if length(self.K) == 1
    self.K = [self.K self.K];
end

%% Capture current figure in high resolution
if ~strcmp(self.figmode,'lazyupdate');
    tempfile = 'lead_temp_screendump.png';
    self.source_fig = gcf;
    current_paperpositionmode = get(self.source_fig,'PaperPositionMode');
    current_inverthardcopy = get(self.source_fig,'InvertHardcopy');
    set(self.source_fig,'PaperPositionMode','auto');
    set(self.source_fig,'InvertHardcopy','off');
    print(self.source_fig,['-r',num2str(1.5*screen_DPI*self.K(1))], '-dpng', tempfile); % change 1.5 to e.g. 2 to get even higher resolution image out.
    set(self.source_fig,'InvertHardcopy',current_inverthardcopy);
    set(self.source_fig,'PaperPositionMode',current_paperpositionmode);
    self.raw_hires = imread(tempfile);
    delete(tempfile);
end
%% Start filtering to remove aliasing
w = warning;
warning off;
if strcmp(self.aamethod,'standard') || strcmp(self.aamethod,'noshrink')
    % Subsample hires figure image with standard anti-aliasing using a
    % butterworth filter
    kk = lpfilter(self.K(2)*3,self.K(2)*0.9,2);
    mm = myconv(ones(size(self.raw_hires(:,:,1))),kk,'same');
    a1 = max(min(myconv(single(self.raw_hires(:,:,1))/(256),kk,'same'),1),0)./mm;
    a2 = max(min(myconv(single(self.raw_hires(:,:,2))/(256),kk,'same'),1),0)./mm;
    a3 = max(min(myconv(single(self.raw_hires(:,:,3))/(256),kk,'same'),1),0)./mm;
    if strcmp(self.aamethod,'standard')
        if abs(1-self.K(2)) > 0.001
            raw_lowres = double(cat(3,a1(2:self.K(2):end,2:self.K(2):end),a2(2:self.K(2):end,2:self.K(2):end),a3(2:self.K(2):end,2:self.K(2):end)));
        else
            raw_lowres = self.raw_hires;
        end
    else
        raw_lowres = double(cat(3,a1,a2,a3));
    end
elseif strcmp(self.aamethod,'imresize')
    % This is probably the fastest method available at this moment...
    raw_lowres = single(imresize(self.raw_hires,1/self.K(2),'bilinear'))/256;
end
warning(w);

%% Place the anti-aliased image in some image on the screen ...
if strcmp(self.figmode,'figure');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    oldpos = get(gcf,'Position');
    self.myaa_figure = figure('Name','Export','Visible','off');
    fig = self.myaa_figure;
    set(fig,'Menubar','none');
    set(fig,'Resize','off');
    sz = size(raw_lowres);
    set(fig,'Units','pixels');
    pos = [oldpos(1:2) sz(2:-1:1)];
    set(fig,'Position',pos);
    ax = axes;
    hi = image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
elseif strcmp(self.figmode,'publish');
    % Create a new figure at the same place as the previous
    % The content of this new image is just a bitmap...
    self.myaa_figure = figure('Name','Export','Visible','off');
    fig = self.myaa_figure;
    current_units = get(self.source_fig,'Units');
    set(self.source_fig,'Units','pixels');
    pos = get(self.source_fig,'Position');
    set(self.source_fig,'Units',current_units);
    set(fig,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    ax = axes;
    hi=image(raw_lowres);
    set(ax,'Units','normalized');
    set(ax,'Position',[0 0 1 1]);
    axis off;
    close(self.source_fig);
elseif strcmp(self.figmode,'update');
    fig = self.myaa_figure;
    figure(fig);
    clf;
    set(fig,'Menubar','none');
    set(fig,'Resize','off');
    sz = size(raw_lowres);
    set(fig,'Units','pixels');
    pos = get(fig,'Position');
    pos(3:4) = sz(2:-1:1);
    set(fig,'Position',pos);
    ax = axes;
    hi=image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
elseif strcmp(self.figmode,'lazyupdate');
    clf;
    fig = self.myaa_figure;
    sz = size(raw_lowres);
    pos = get(fig,'Position');
    pos(3:4) = sz(2:-1:1);
    set(fig,'Position',pos);
    ax = axes;
    hi=image(raw_lowres);
    set(ax,'Units','pixels');
    set(ax,'Position',[1 1 sz(2) sz(1)]);
    axis off;
end

%% Store current state

set(gcf,'userdata',self);
set(gcf,'KeyPressFcn',@keypress);
set(gcf,'Interruptible','off');

%% Avoid unnecessary console output
if nargout == 1
    varargout(1) = {fig};
elseif nargout == 2
    varargout(1) = {fig};
    varargout(2) = {get(hi, 'CData')};
    close(self.myaa_figure);
end


%% A simple lowpass filter kernel (Butterworth).
% sz is the size of the filter
% subsmp is the downsampling factor to be used later
% n is the degree of the butterworth filter
function kk = lpfilter(sz, subsmp, n)
sz = 2*floor(sz/2)+1; % make sure the size of the filter is odd
cut_frequency = 0.5 / subsmp;
range = (-(sz-1)/2:(sz-1)/2)/(sz-1);
[ii,jj] = ndgrid(range,range);
rr = sqrt(ii.^2+jj.^2);
kk = ifftshift(1./(1+(rr./cut_frequency).^(2*n)));
kk = fftshift(real(ifft2(kk)));
kk = kk./sum(kk(:));


function keypress(src,evnt)
if isempty(evnt.Character)
    return
end

recognized = 0;
self = get(gcf,'userdata');

if evnt.Character == '+'
    self.K(2) = max(self.K(2).*0.5^(1/2),1);
    recognized = 1;
    set(gcf,'userdata',self);
    myaa('lazyupdate');
elseif evnt.Character == '-'
    self.K(2) = min(self.K(2).*2^(1/2),16);
    recognized = 1;
    set(gcf,'userdata',self);
    myaa('lazyupdate');
elseif evnt.Character == ' ' || evnt.Character == 'r' || evnt.Character == 'R'
    set(gcf,'userdata',self);
    myaa('update');
elseif evnt.Character == 'q'
    close(gcf);
elseif find('123456789' == evnt.Character)
    self.K = [str2double(evnt.Character) str2double(evnt.Character)];
    set(gcf,'userdata',self);
    myaa('update');
end


function ea_distogrotate
uibjs=getappdata(gcf,'uibjs');
set(uibjs.rotate3dtog,'State','off');


function ea_distogzoomin
uibjs=getappdata(gcf,'uibjs');
set(uibjs.magnifyplus,'State','off');


function ea_distogzoomout
uibjs=getappdata(gcf,'uibjs');
set(uibjs.magnifyminus,'State','off');


function ea_distogpan
uibjs=getappdata(gcf,'uibjs');
set(uibjs.handtog,'State','off');


function ea_rotate(h,~,cmd)

h = rotate3d;
h.RotateStyle = 'orbit';
h.Enable = cmd;
if strcmp(cmd,'on')
    ea_distogzoomin;
    ea_distogzoomout;
    ea_distogpan;
end

function ea_pan(h,~,cmd)

pan(cmd);
if strcmp(cmd,'on')
    ea_distogzoomin;
    ea_distogzoomout;
    ea_distogrotate;
end


function ea_zoomin(h,~,cmd)

if strcmp(cmd,'on')
    ea_distogpan;
    ea_distogrotate;
    ea_distogzoomout;
end
h=zoom;
h.Enable=cmd;
h.Motion='both';
h.Direction='in';


function ea_zoomout(h,~,cmd)

if strcmp(cmd,'on')
    ea_distogpan;
    ea_distogrotate;
    ea_distogzoomin;
end
h=zoom;
h.Enable=cmd;
h.Motion='both';
h.Direction='out';

function [handle,string] = ea_setmertag(tag,keymer,trajectory)
% tag.string; tag.visible; tag.color;
d = 3.2;
pos = [trajectory(1,1)/abs(trajectory(1,1))*d+trajectory(1,1),trajectory(1,2:3)];
% string = sprintf('%s%s Depth: %smm',upper(keymer(4)),keymer(5:strfind(keymer,'_')-1),tag.depth);
string = sprintf('%s%s: %smm',upper(keymer(4)),keymer(5:strfind(keymer,'_')-1),tag.depth);
handle = text(pos(1),pos(2),pos(3),string,'Color','w','HorizontalAlignment','center','Visible',tag.visible);


function ea_updatemercontrol(keymer,handles,mermarkers,resultfig,options)
n=length(mermarkers);
markerstring.right = get(handles.popupmermarkers_right,'String');
markerstring.left = get(handles.popupmermarkers_left,'String');

if isempty(markerstring.right)
    markerstring.right = {'none selected...'};
end

if isempty(markerstring.left)
    markerstring.left = {'none selected...'};
end

if strcmp(keymer(strfind(keymer,'_')+1:end),'right')
    % side = 1;
    markerstring.right{end+1} = sprintf('%0.0f. %s',n,mermarkers(n).tag.string);
elseif strcmp(keymer(strfind(keymer,'_')+1:end),'left')
    % side = 2;
    markerstring.left{end+1} = sprintf('%0.0f. %s',n,mermarkers(n).tag.string);
end

setappdata(resultfig,'markerstring',markerstring)
set(handles.popupmermarkers_right,'Visible','off','String',markerstring.right,'Value',1)
set(handles.popupmermarkers_left,'Visible','off','String',markerstring.left,'Value',1)
% set(handles.togglemarkertags,'Visible','on','Value',1)

function outputtrajectory = ea_getmertrajectory(trajectory,dist,length,n)
if size(trajectory,1)<2
    error('Must input a vector')
end
dxyz = sqrt((diff(trajectory(1:2,1))^2)+(diff(trajectory(1:2,2))^2)+diff(trajectory(1:2,3))^2);
slope = mean(diff(trajectory))/dxyz;
startpoint = trajectory(1,:)+slope.*dist;

outputtrajectory(:,1) = linspace(startpoint(1,1),startpoint(1,1)+slope(1)*length,n);
outputtrajectory(:,2) = linspace(startpoint(1,2),startpoint(1,2)+slope(2)*length,n);
outputtrajectory(:,3) = linspace(startpoint(1,3),startpoint(1,3)+slope(3)*length,n);


function ea_updatemertrajectory(handles,trajectory,dist,tag) 
resultfig=getappdata(handles.mercontrolfig,'resultfig');
% Update position in resultfig
% XData = get(getappdata(resultfig,tag),'XData');
% YData = get(getappdata(resultfig,tag),'YData');
% ZData = get(getappdata(resultfig,tag),'ZData');
h = getfield(getappdata(resultfig,'merhandles'),tag);
set(h,'XData',trajectory(:,1)')
set(h,'YData',trajectory(:,2)')
set(h,'ZData',trajectory(:,3)')
setappdata(resultfig,tag,h)
set(handles.(['key',tag]),'Value',1)
setappdata(resultfig,'keymer',['key',tag])
newdiststr = num2str(str2double(get(handles.(['pos',tag]),'String'))+dist);
set(handles.(['pos',tag]),'String',newdiststr)