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

if ~isfield(options, 'leadprod')
    options.leadprod = '';
end

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

titlePrefix = erase(options.patientname, 'sub-');

resultfig=figure('name', [titlePrefix,': Electrode-Scene'],...
    'color', 'k', 'numbertitle', 'off',...
    'CloseRequestFcn', @closesatellites, 'visible', options.d3.verbose,...
    'KeyPressFcn', @ea_keypress, 'KeyReleaseFcn', @ea_keyrelease);

setappdata(resultfig,'options',options);

ea_bind_dragndrop(resultfig, ...
    @(obj,evt) DropFcn(obj,evt,resultfig), ...
    @(obj,evt) DropFcn(obj,evt,resultfig));

set(resultfig,'toolbar','none');
ssz=get(0,'Screensize');
ssz(1:2)=ssz(1:2)+50;
ssz(3:4)=ssz(3:4)-200;
set(resultfig, 'Position', ssz); % Maximize figure.
options.d3.exportBB=0;
% initialize some ui elements
ht=uitoolbar(resultfig);
setappdata(resultfig,'ht',ht);

% add custom rotator:
uibjs.rotate3dtog=uitoggletool(ht, 'CData', ea_get_icn('rotate'),...
    'TooltipString', 'Rotate - Pan - Zoom', 'OnCallback', {@ea_rotate,'on'},...
    'OffCallback', {@ea_rotate,'off'}, 'State', 'on');
uibjs.slide3dtog=uitoggletool(ht, 'CData', ea_get_icn('quiver'),...
    'TooltipString', 'Slide Slices', 'OnCallback', {@ea_slideslices,'on'},...
    'OffCallback', {@ea_slideslices,'off'}, 'State', 'off');
% uibjs.magnifyplus=uitoggletool(ht,'CData',ea_get_icn('magnplus'),...
%     'TooltipString', 'Zoom In', 'OnCallback', {@ea_zoomin,'on'},...
%     'OffCallback', {@ea_zoomin,'off'}, 'State', 'off');
% uibjs.magnifyminus=uitoggletool(ht, 'CData', ea_get_icn('magnminus'),...
%     'TooltipString', 'Zoom Out', 'OnCallback', {@ea_zoomout,'on'},...
%     'OffCallback', {@ea_zoomout,'off'}, 'State', 'off');
% uibjs.handtog=uitoggletool(ht, 'CData', ea_get_icn('hand'),...
%     'TooltipString', 'Pan Scene', 'OnCallback', {@ea_pan,'on'},...
%     'OffCallback', {@ea_pan,'off'}, 'State', 'off');
setappdata(resultfig,'uibjs',uibjs);

% Initialize Sliceview-Button
slicebutton=uipushtool(ht,'CData',ea_get_icn('slices'),...
    'TooltipString','Slice Control Figure',...
    'ClickedCallback',{@opensliceviewer,resultfig,options});

mh = uimenu(resultfig,'Label','Add Objects');
fh1 = uimenu(mh,'Label','Open Tract',...
    'Callback',{@(src, evt) ea_addobj(resultfig,'tract',options)});
fh2 = uimenu(mh,'Label','Open ROI',...
    'Callback',{@(src, evt) ea_addobj(resultfig,'roi',options)});
fh3 = uimenu(mh,'Label','Show tracts weighted by ROI',...
    'Callback',{@(src, evt) ea_addobj(resultfig,'tractmap',options)});
fh3 = uimenu(mh,'Label','Show fiber activation result from OSS-DBS',...
    'Callback',{@(src, evt) ea_addobj(resultfig,'fiberactivation',options)});

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
prefs=ea_prefs;
% colormap('gray')

%% Patient specific part (skipped if no patient is selected or no reco available):
if ~strcmp(options.patientname,'No Patient Selected') % if not initialize empty viewer
    if nargin>1 || isfield(options.subj, 'recon') && isfile(options.subj.recon.recon)
        if nargin>1
            multiplemode=1;

            % mer development
            % if isstruct(varargin{2})
            %     elstruct=varargin{2}.elstruct;
            %     merstruct=varargin{2}.merstruct;
            % else
            elstruct=varargin{2};
            % end

            if options.d3.mirrorsides
                elstruct=ea_mirrorsides(elstruct);
                try options.d3.isomatrix=ea_mirrorsides(options.d3.isomatrix); end
            end
        else
            multiplemode=0;
            [coords_mm,trajectory,markers]=ea_load_reconstruction(options);
            elstruct(1).coords_mm=coords_mm;
            elstruct(1).trajectory=trajectory;
            elstruct(1).name=options.patientname;
            elstruct(1).markers=markers;
            clear coords_mm trajectory
        end

        elSide = cell(1, length(elstruct));
        for pt=1:length(elstruct)
            if exist('el_render','var')
                [el_render,el_label,elSide{pt}]=ea_renderelstruct(options,resultfig,elstruct,pt,el_render,el_label);
            else
                [el_render,el_label,elSide{pt}]=ea_renderelstruct(options,resultfig,elstruct,pt);
            end

            if ~multiplemode
                side=options.sides(end);
                d=load(options.subj.recon.recon);
                plans=d.reco.electrode(side+1:end);
                if ~isempty(plans)
                    if isfield(plans,'plan')
                        for plan=1:length(plans)
                            pobj=ea_load_electrode(options.subj.recon.recon, side+plan);
                            ea_add_trajectory([],[],options,pobj,side+plan);
                        end
                    end
                end
                eltext=getappdata(resultfig,'eltext');

                eltexttoggle=uitoggletool(ht, 'CData', ea_get_icn('electrode_segment'),...
                    'TooltipString', 'Contact Labels',...
                    'OnCallback', {@objvisible,eltext},...
                    'OffCallback', {@objinvisible,eltext}, 'State','off');


            end
            if options.d3.elrendering==1 && options.d3.exportBB % export vizstruct for lateron export to JSON file / Brainbrowser.
                % this part for brainbrowser support.
                vizstruct=struct('faces',[],'vertices',[],'colors',[]);

                cnt=1;
                for iside=1:length(options.sides)
                    side=options.sides(iside);
                    extract=1:length(el_render(side).elpatch);
                    for ex=extract
                        tp=el_render(side).elpatch(ex);

                        try % works only in ML 2015
                            tr=triangulation(get(el_render(side).elpatch(ex),'Faces'),...
                                get(el_render(side).elpatch(ex),'Vertices'));
                            vizstruct(cnt).normals = vertexNormal(tr);
                        catch % workaround for older versions..
                            vizstruct(cnt).normals=get(tp,'VertexNormals')*-1;
                        end

                        vizstruct(cnt).faces=get(tp,'Faces');
                        vizstruct(cnt).vertices=get(tp,'Vertices');
                        scolor=get(el_render(side).elpatch(ex),'FaceVertexCData');
                        vizstruct(cnt).colors=scolor;
                        %vizstruct(cnt).colors=repmat([squeeze(scolor(1,1,:))',0.7],length(vizstruct(cnt).faces),1);
                        vizstruct(cnt).name='';
                        cnt=cnt+1;
                    end
                end
            end

            % show microelectrode recording data
            if exist('merstruct','var')
                try
                    [mer(pt).render,merlabel(:,pt)]=ea_showmer(resultfig,merstruct(pt),pt,options);
                catch
                    ea_error(['Couldn''''t visualize electrode from patient ',num2str(pt),'.']);
                end
            end
        end

        setappdata(resultfig,'elstruct',elstruct);
        setappdata(resultfig,'el_render',el_render);
        % add handles to buttons. Can't be combined with the above loop since all
        % handles need to be set for the buttons to work properly (if alt is
        % pressed, all electrodes are made visible/invisible).
        drawnow

        if strcmp(options.leadprod,'group')
            elstructGroupID = [elstruct.group];
            sideNum = cellfun(@numel, elSide);
            elrenderGroupID = repelem(elstructGroupID, sideNum);
            for g = unique(elstructGroupID)
                el_renderID = elrenderGroupID == g;
                eleGroupToggle(g) = uitoggletool(ht, 'CData', ea_get_icn('electrode_group'),...
                    'TooltipString', ['Electrode Group ', num2str(g)],...
                    'Tag', ['Group: ', num2str(g)],...
                    'OnCallback', {@eleGroupVisible,el_render(el_renderID)},...
                    'OffCallback', {@eleGroupInvisible,el_render(el_renderID)}, 'State','on');
            end
            setappdata(resultfig,'eleGroupToggle',eleGroupToggle);

            % add sweetspot explorer button.
            sweetspotadd = uipushtool(ht, 'CData', ea_get_icn('sweetspot_add'),...
                'TooltipString', ['Add sweetspot analysis'],...
                'Tag', ['Add sweetspot analysis'],...
                'ClickedCallback', {@ea_add_sweetspot,ea_getGroupAnalysisFile(options.groupdir),resultfig});

            di=dir([options.root,options.patientname,filesep,'sweetspots',filesep,'*.sweetspot']);
            for d=1:length(di)
                uipushtool(ht, 'CData', ea_get_icn('sweetspot'),...
                    'TooltipString', ['Explore sweetspot analysis ',ea_stripext(di(d).name)],...
                    'Tag', ['Explore sweetspot analysis ',ea_stripext(di(d).name)],...
                    'ClickedCallback', {@ea_add_sweetspot,fullfile(options.groupdir,'sweetspots',di(d).name),resultfig});
            end


            % add discriminative fiber explorer button.
            discfiberadd = uipushtool(ht, 'CData', ea_get_icn('discfiber_add'),...
                'TooltipString', ['Add fiber filtering analysis'],...
                'Tag', ['Add fiber filtering analysis'],...
                'ClickedCallback', {@ea_add_discfiber,ea_getGroupAnalysisFile(options.groupdir),resultfig});

            di=dir([options.root,options.patientname,filesep,'fiberfiltering',filesep,'*.fibfilt']);
            for d=1:length(di)
                uipushtool(ht, 'CData', ea_get_icn('discfiber'),...
                    'TooltipString', ['Explore fiber filtering analysis ',ea_stripext(di(d).name)],...
                    'Tag', ['Explore fiber filtering analysis ',ea_stripext(di(d).name)],...
                    'ClickedCallback', {@ea_add_discfiber,fullfile(options.groupdir,'fiberfiltering',di(d).name),resultfig});
            end

            % add networkmapping explorer button.
            netmapadd = uipushtool(ht, 'CData', ea_get_icn('networkmapping_add'),...
                'TooltipString', ['Add DBS network mapping analysis'],...
                'Tag', ['Add DBS network mapping analysis'],...
                'ClickedCallback', {@ea_add_networkmapping,ea_getGroupAnalysisFile(options.groupdir),resultfig});

            di=dir([options.root,options.patientname,filesep,'networkmapping',filesep,'*.netmap']);
            for d=1:length(di)
                uipushtool(ht, 'CData', ea_get_icn('networkmapping'),...
                    'TooltipString', ['Explore DBS network mapping analysis ',ea_stripext(di(d).name)],...
                    'Tag', ['Explore DBS network mapping analysis ',ea_stripext(di(d).name)],...
                    'ClickedCallback', {@ea_add_networkmapping,fullfile(options.groupdir,'networkmapping',di(d).name),resultfig});
            end


            % Move the group toggles and app toggles forward
            isEleToggle = arrayfun(@(obj) ~isempty(regexp(obj.Tag, '^Group: \d+,', 'once')), allchild(ht));
            eleToggleInd = find(isEleToggle);
            isEleGroupToggle = arrayfun(@(obj) ~isempty(regexp(obj.Tag, '^Group: \d+$', 'once')), allchild(ht));
            eleGroupToggleInd = find(isEleGroupToggle);
            otherToggleInd = (find(isEleToggle,1,'last')+1:numel(ht.Children))';
            appToggleInd = (1:find(isEleGroupToggle,1)-1)';
            ht.Children=ht.Children([eleToggleInd;eleGroupToggleInd;appToggleInd;otherToggleInd]);
        end

        try
            set(el_label,'Visible','off');
            ellabeltog = uitoggletool(ht, 'CData', ea_get_icn('labels'),...
                'TooltipString', 'Electrode Labels',...
                'OnCallback', {@objvisible,el_label},...
                'OffCallback', {@objinvisible,el_label}, 'State','off');

            % Move eleLabel toggle to front
            if strcmp(options.leadprod,'dbs')
                eleToggleTagPattern = '^Patient: ';
            elseif strcmp(options.leadprod,'group')
                eleToggleTagPattern = '^Group: \d+';
            end
            isEleToggle = arrayfun(@(obj) ~isempty(regexp(obj.Tag, eleToggleTagPattern, 'once')), allchild(ht));
            eleToggleInd = find(isEleToggle);
            otherToggleInd = (find(isEleToggle,1,'last')+1:numel(ht.Children))';
            ht.Children=ht.Children([eleToggleInd;1;otherToggleInd]);
        end

        cnt=1;
        for pt=1:length(elstruct)
            try
                if multiplemode
                    caption{1}=[elstruct(pt).name,'_Left'];
                    caption{2}=[elstruct(pt).name,'_Right'];
                else
                    caption{1}='Electrode_Left';
                    caption{2}='Electrode_Right';
                end
                %eltog(cnt)=uitoggletool(ht,'CData',ea_get_icn('electrode'),'TooltipString',caption{1},'OnCallback',{@elvisible,el_render,pt,2,'on',options},'OffCallback',{@elvisible,el_render,pt,2,'off',options},'State','on');
                %eltog(cnt+1)=uitoggletool(ht,'CData',ea_get_icn('electrode'),'TooltipString',caption{2},'OnCallback',{@elvisible,el_render,pt,1,'on',options},'OffCallback',{@elvisible,el_render,pt,1,'off',options},'State','on');
                if isfield(options,'uipatdirs')
                    if exist([options.uipatdirs{pt} '/cortex/CortElecs.mat'],'file')
                        vars = whos('-file',[options.uipatdirs{pt} '/cortex/CortElecs.mat']);
                        CortElecs = load([options.uipatdirs{pt} '/cortex/CortElecs.mat']);
                        if ismember('Left',{vars.name})
                            hold on; plot3(CortElecs.Left(:,1),CortElecs.Left(:,2),CortElecs.Left(:,3),'.','color','r','markersize',10)
                            ctxeltog(cnt)=uitoggletool(ht,'CData',ea_get_icn('cortical_strip'),'TooltipString',['Cortical_' caption{1}],'OnCallback',{@ctxelvisible,el_renderstrip,pt,2,'on',options},'OffCallback',{@elvisible,el_render,pt,2,'off',options},'State','on');
                        end
                        if ismember('Right',{vars.name})
                            ctxeltog(cnt)=uitoggletool(ht,'CData',ea_get_icn('cortical_strip'),'TooltipString',['Cortical_' caption{1}],'OnCallback',{@ctxelvisible,el_renderstrip,pt,2,'on',options},'OffCallback',{@elvisible,el_render,pt,2,'off',options},'State','on');
                        end
                    end
                end
                cnt=cnt+2;
            end
        end

        clear cnt

        % Initialize Stimulation-Button
        if ~strcmp(options.leadprod, 'group')
            eladdTraj = uipushtool(ht,'CData',ea_get_icn('addelectrode'),...
                'TooltipString','Add Trajectory...','ClickedCallback',{@ea_add_trajectory,options});
            stimbutton = uipushtool(ht,'CData',ea_get_icn('stimulation'),...
                'TooltipString','Stimulation Control Figure',...
                'ClickedCallback',{@openstimviewer,elstruct,resultfig,options});
        end

    else
        options.writeoutstats=0; % if no electrodes are there, stats can't be written.
        elstruct=struct;
    end

else
    options.writeoutstats=0; % if no electrodes are there, stats can't be written.
    elstruct=struct;
end

% Initialize MER-Button
if ~strcmp(options.leadprod, 'group')
    merbutton=uipushtool(ht,'CData',ea_get_icn('mer'),...
        'TooltipString','MER Control Figure',...
        'ClickedCallback',{@ea_openmerviewer,resultfig,options});
end

% Initialize Convis-Button
if ~strcmp(options.leadprod,'group') 
convisbutton=uipushtool(ht,'CData',ea_get_icn('connectome'),...
    'TooltipString','Connectivity Visualization',...
    'ClickedCallback',{@openconnectomeviewer,resultfig,options});
end

% Initialize FS Cortex-Button
corticalbutton=uipushtool(ht,'CData',ea_get_icn('cortex'),...
    'TooltipString','Cortical Reconstruction Visualization',...
    'ClickedCallback',{@opencortexviewer,resultfig,options});

% Initialize Cortical Strip-Button
% cortelsbutton=uipushtool(ht,'CData',ea_get_icn('cortical_strip'),...
%     'TooltipString','Cortical Reconstruction Visualization',...
%     'ClickedCallback',{@opencortelsviewer,resultfig,options});

% set defaultview
prefs = ea_prefs;
v = prefs.machine.view;
ea_view(v);

% Show atlas data
if options.d3.writeatlases && ~strcmp(options.atlasset, 'Use none')
    atlases = ea_showatlas(resultfig,elstruct,options);

    if ~strcmp(options.d3.verbose,'off') && ~atlases.discfibersonly
        ea_openatlascontrol([],[],atlases,resultfig,options);
    end
else
    colormap(gray)
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

if isfield(options.d3,'expdf')
    if options.d3.expdf
        exportDir = [options.subj.exportDir, filesep, 'pdf'];
        ea_mkdir(exportDir);
        fig2pdf3d(gca,[exportDir,filesep,'Lead-DBS_Electrode_Localization'],options);
        close(resultfig);
        return
    end
end
%% End of patient's part.

% Initialize a draggable lightbulb
set(0,'CurrentFigure',resultfig);
hold on
ea_show_light(resultfig,1);
% set(lightbulb, 'Visible', 'off');

lightbulbbutton=uipushtool(ht,'CData',ea_get_icn('lightbulb'),...
    'TooltipString','Set Lighting',...
    'ClickedCallback',{@ea_launch_setlighting,resultfig});
% clightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('clightbulb'),...
%     'TooltipString','Ceiling Lightbulb',...
%     'OnCallback',{@objvisible,getappdata(resultfig,'ceiling_lamp')},...
%     'OffCallback',{@objinvisible,getappdata(resultfig,'ceiling_lamp')},'State','on');
% llightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('llightbulb'),...
%     'TooltipString','Left Lightbulb',...
%     'OnCallback',{@objvisible,getappdata(resultfig,'left_lamp')},...
%     'OffCallback',{@objinvisible,getappdata(resultfig,'left_lamp')},'State','on');
% rlightbulbbutton=uitoggletool(ht,'CData',ea_get_icn('rlightbulb'),...
%     'TooltipString','Right Lightbulb',...
%     'OnCallback',{@objvisible,getappdata(resultfig,'right_lamp')},...
%     'OffCallback',{@objinvisible,getappdata(resultfig,'right_lamp')},'State','on');

if options.prefs.env.dev
    setBlackBackgroundButton = uipushtool(ht,'CData',ea_get_icn('BGB'),...
        'TooltipString','Set background to Black',...
        'ClickedCallback',@(~, ~) ea_setElvisBlackBackground(resultfig));
    setWhiteBackgroundButton = uipushtool(ht,'CData',ea_get_icn('BGW'),...
        'TooltipString','Set background to White',...
        'ClickedCallback',@(~, ~) ea_setElvisWhiteBackground(resultfig));
    setTransparentBackgroundButton = uipushtool(ht,'CData',ea_get_icn('BGT'),...
        'TooltipString','Set background to Transparent',...
        'ClickedCallback',@(~, ~) ea_setElvisTransparentBackground(resultfig));
end

% Initialize HD-Export button
dumpscreenshotbutton=uipushtool(ht,'CData',ea_get_icn('dump'),...
    'TooltipString','Dump Screenshot','ClickedCallback',{@dump_screenshot,resultfig,options});

hdsavebutton=uipushtool(ht,'CData',ea_get_icn('save'),...
    'TooltipString','Save Scene','ClickedCallback',@export_hd);
dofsavebutton=uipushtool(ht,'CData',ea_get_icn('save_depth'),...
    'TooltipString','Save Scene with depth of field',...
    'ClickedCallback',{@ea_export_depth_of_field,resultfig});

% Initialize Video-Export button
videoexportbutton=uipushtool(ht,'CData',ea_get_icn('video'),...
    'TooltipString','Save video','ClickedCallback',{@export_video,options});

% Init hard_electrode_view button
if isfield(options,'modality') && options.modality==2
    electrodesegmentbutton=uitoggletool(ht,'CData',ea_get_icn('electrode_segment'),...
        'TooltipString','Auto-Segment electrode from postoperative acquisition',...
        'OnCallback',{@ea_segment_electrode,options,resultfig,'on'},...
        'OffCallback',{@ea_segment_electrode,options,resultfig,'off'},'State','off');
end
% Initialize Export to Lead-Server button

% lsbutton=uipushtool(ht,'CData',ea_get_icn('server'),...
%     'TooltipString','Export to Server',...
%     'ClickedCallback',{@ea_export_server,options});

% add default view buttons
uipushtool(ht, 'CData',ea_get_icn('defaultviewsave'),...
    'TooltipString', 'Save current view as default.',...
    'ClickedCallback',@save_currentview_callback);
uipushtool(ht, 'CData',ea_get_icn('defaultviewset'),...
    'TooltipString', 'Display default view',...
    'ClickedCallback',@set_defaultview_callback);

% Reorder toggles to move eleToggle to the end
if strcmp(options.leadprod,'group')
    isEleToggle = arrayfun(@(obj) ~isempty(regexp(obj.Tag, '^Group: \d+', 'once')), allchild(ht));
    eleToggleInd = find(isEleToggle);
    eleLabelToggleInd = eleToggleInd(end) + 1;
    otherToggleInd = setdiff((1:numel(ht.Children))', [eleToggleInd;eleLabelToggleInd]);
    ht.Children=ht.Children([eleToggleInd;eleLabelToggleInd;otherToggleInd]);
end

hold off

set(0,'CurrentFigure',resultfig);

set(resultfig,'Renderer','OpenGL')
axis off
set(resultfig,'color','k');

axis vis3d
axis equal
set(resultfig,'Name',figtitle);
set(0,'CurrentFigure',resultfig);

try
    set(gca,'cameraviewanglemode','manual');
end
set(gca,'clipping','off');

%set(resultfig,'visible','on');
if ~strcmp(options.d3.verbose,'off')
    opensliceviewer([],[],resultfig,options);
end

try
    togglestates = prefs.machine.togglestates;
    ea_defaultview(v,togglestates);
end

if options.d3.elrendering==1 && options.d3.exportBB % export vizstruct for lateron export to JSON file / Brainbrowser.
    try
        % store json in figure file
        bbstruct=ea_viz2brainbrowser(vizstruct);
        setappdata(resultfig,'bbstruct',bbstruct);
    end
    if options.prefs.ls.autosave
        ea_export_server([],[],options);
    end
end

setappdata(resultfig, 'options', options);
setappdata(resultfig,'elstruct',elstruct);
ea_figmenu(resultfig,'add');

% set mouse camera opts
ax = findobj(resultfig.Children,'Type','axes');
set(findobj(ax.Children,'Type','surface'),'HitTest','off');
ea_mouse_camera(resultfig);


function ea_launch_setlighting(~,~,resultfig)
ea_set_lighting(resultfig);

% --- Drag and drop callback to load patdir.
function DropFcn(~, event, resultfig)

switch event.DropType
    case 'file'
        objects = event.Data;
    case 'string'
        objects = {event.Data};
end

nonexist = cellfun(@(x) ~exist(x, 'file'), objects);
if any(nonexist)
    fprintf('\nExcluded non-existent/invalid folder:\n');
    cellfun(@disp, objects(nonexist));
    fprintf('\n');
    objects(nonexist) = [];
end

options = getappdata(resultfig,'options');

if ~isempty(objects)
    ea_addobj(resultfig, objects, options);
end


function opensliceviewer(hobj,ev,resultfig,options)
awin=ea_anatomycontrol(resultfig,options);

set(awin,'visible',options.d3.verbose);
setappdata(resultfig,'awin',awin);


function openconnectomeviewer(hobj,ev,resultfig,options)
conwin=ea_convis(gcf,options);
setappdata(resultfig,'conwin',conwin);


function openstimviewer(hobj,ev,elstruct,resultfig,options)
stimwin=ea_stimparams(elstruct,gcf,options);
setappdata(resultfig,'stimwin',stimwin);
%try WinOnTop(stimwin,true); end


function opencortexviewer(hobj,ev,resultfig,options)
cortex=ea_showcortex(resultfig,options);
setappdata(resultfig,'cortex',cortex);
% reload slice viewer to update opacity control
awin=ea_anatomycontrol(resultfig,options);
setappdata(resultfig,'awin',awin);
try WinOnTop(awin,true); end


function closesatellites(src,evnt)
% Close all valid satellite windows
structfun(@(f) isa(f, 'matlab.ui.Figure') && isvalid(f) && close(f), getappdata(gcf));
delete(gcf)


% default view buttons callback
function save_currentview_callback(source,eventdata)
% call ea_defaultview so current view is saved
ea_defaultview()


function set_defaultview_callback(source,eventdata)
% get stored default view preferences and call ea_defaultview
prefs = ea_prefs;
v = prefs.machine.view;
togglestates = prefs.machine.togglestates;
ea_defaultview_transition(v,togglestates);
ea_defaultview(v,togglestates);


function ea_setElvisBlackBackground(resultfig)
cmap = gray;
set(resultfig, 'Color', 'k', 'Colormap', cmap);


function ea_setElvisWhiteBackground(resultfig)
cmap = gray;
% Get volume data
V = getappdata(resultfig, 'V');
if isa(V{1}, 'nifti')
    V = V{1}.dat; % Memory mapped nifti struct
else
    V = V{1}.img; % Standard nifti struct
end

% Take the middle z slice
zslice = V(:,:,round(size(V,3)/2));

% Check if background (1st voxel) is dark or bright
if zslice(1,1) < mean(zslice(:))
    % Flip black to white in colormap in case background is dark
    cmap(1,:) = [1 1 1];
end

set(resultfig, 'Color', 'w', 'Colormap', cmap);


function ea_setElvisTransparentBackground(resultfig)
cmap = gray;
% Get volume data
V = getappdata(resultfig, 'V');
if isa(V{1}, 'nifti')
    V = V{1}.dat; % Memory mapped nifti struct
else
    V = V{1}.img; % Standard nifti struct
end

% Take the middle z slice
zslice = V(:,:,round(size(V,3)/2));

% Check if background (1st voxel) is dark or bright
if zslice(1,1) < mean(zslice(:))
    % Flip black to white in colormap in case background is dark
    cmap(1,:) = [1 1 1];
end

set(resultfig, 'Color', 'none', 'Colormap', cmap);


function export_video(hobj,ev,options)
% Set up recording parameters (optional), and record
[FileName,PathName] = uiputfile('LEAD_Scene.mp4','Save file name for video');
ea_CaptureFigVid(options.prefs.video.path, [PathName,FileName],options.prefs.video.opts);


function export_hd(hobj,ev)
set(gca, 'Color', 'none');
set(gcf,'color','none');

[FileName,PathName] = uiputfile('LEAD_Scene.png','Save file name');
if FileName
    % set(gcf, 'Color', [1,1,1]);
    % [~, cdata] = ea_myaa([4, 2]);
    %
    % imwrite(cdata, [PathName,FileName], 'png');
    ea_screenshot([PathName,FileName],'myaa');
end


function dump_screenshot(hobj,ev,resultfig,options)

set(0,'CurrentFigure',resultfig);
if ~exist([options.root,options.patientname,filesep,'export',filesep,'views'],'dir')
    mkdir([options.root,options.patientname,filesep,'export',filesep,'views']);
end

views=dir([options.root,options.patientname,filesep,'export',filesep,'views',filesep,'view_u*.png']);
next=1;
while ismember(['view_u',sprintf('%03.0f',next),'.png'],{views.name})
    next=next+1;
end

ea_screenshot([options.root,options.patientname,filesep,'export',filesep,'views',filesep,'view_u',sprintf('%03.0f',next),'.png'],'ld');


function eleGroupVisible(hobj,ev,eleGroup)
for i=1:numel(eleGroup)
    eleGroup(i).toggleH.State='on';
end


function eleGroupInvisible(hobj,ev,eleGroup)
for i=1:numel(eleGroup)
    eleGroup(i).toggleH.State='off';
end


function objvisible(hobj,ev,atls)
set(atls, 'Visible', 'on');


function objinvisible(hobj,ev,atls)
set(atls, 'Visible', 'off');


function ctxelvisible(hobj,ev,atls,pt,side,onoff,options)

if(getappdata(gcf,'altpressed'))

    eltog=getappdata(hobj.Parent.Parent,'eltog');
    set(eltog,'State',onoff);
    for el=1:length(atls)
        for iside=1:length(options.sides)
            side=options.sides(iside);
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


function ea_keypress(resultfig, event)
% this is the main keypress function for the resultfigure. Add event
% listeners here.
if ismember('alt', event.Modifier)
    setappdata(resultfig, 'altpressed', 1);
elseif ismember('shift', event.Modifier)
    setappdata(resultfig, 'shiftpressed', 1);
elseif ismember('command', event.Modifier)
    setappdata(resultfig, 'cmdpressed', 1);
end

% If the MER Control window is open
mercontrolfig = getappdata(resultfig, 'mercontrolfig');
if ~isempty(mercontrolfig) && isvalid(mercontrolfig)

    merstruct = getappdata(mercontrolfig, 'merstruct');

    bChecked = logical([merstruct.Toggles.keycontrol.value]);
    if ~any(bChecked)
        return
    end

    if any(strcmpi(event.Key, {'space','m','l','t','b', 's', 'n'}))
        % Reserved keys:
        % 'space' = Generic; 'm' = MER; 'l' = LFP; 't' = Top; 'b' = Bottom
        % 's' = session; 'n' = notes

        if any(strcmpi(event.Key, {'s', 'n'}))
            % Enter session or notes for the last marker.
            if strcmpi(event.Key, 's')
                merstruct.markers(end).session = char(inputdlg('Enter Session'));
            elseif strcmpi(event.Key, 'n')
                merstruct.markers(end).notes = char(inputdlg('Enter Notes'));
            end
            setappdata(resultfig, 'mermarkers', merstruct.markers);
            return;
        end

        sess_text = '';
        switch lower(event.Key)
            case 'space'
                markertype = MERState.MarkerTypes.Generic;
            case 'm'
                markertype = MERState.MarkerTypes.MER;
                sess_text = char(inputdlg ('Enter Session'));
            case 'l'
                markertype = MERState.MarkerTypes.LFP;
                sess_text = char(inputdlg ('Enter Session'));
            case 't'
                markertype = MERState.MarkerTypes.Top;
            case 'b'
                markertype = MERState.MarkerTypes.Bottom;
        end

        % For each checked box, add a marker.
        merstruct.addMarkersAtTrajs(merstruct.Toggles.keycontrol(bChecked),...
            markertype, sess_text);
        handles = guidata(mercontrolfig);
        ea_resultfig_updatemarkers(handles);
        ea_mercontrol_updatemarkers(handles);

    elseif any(strcmpi(event.Key, {'uparrow','leftarrow','downarrow','rightarrow'}))
        d = 1;  % Default step size
        if ~isempty(event.Modifier)
            if ismember(event.Modifier, 'shift')
                d = 2;  % Large step size
            elseif ismember(event.Modifier, 'alt')
                d = 3;  % Small step size
            end
        end
        if any(strcmpi(event.Key, {'downarrow','rightarrow'}))
            d = -d;
        end
        merstruct.translateToggledTrajectories(d);

        % Update the GUI
        handles = guidata(mercontrolfig);
        ea_resultfig_updatetrajectories(handles);
        ea_mercontrol_updateimplanted(handles);

    end
    % commnd=event.Character;
    % switch lower(commnd)
    %     case ' '
    %     case {'x','a','p','y','l','r'} % view angles.
    %     case {'0','3','4','7'}
    %     case {'c','v','b','n'}
    %     otherwise % arrow keys, plus, minus
    % end
end


function ea_keyrelease(resultfig,event)
setappdata(resultfig,'altpressed',0);
setappdata(resultfig,'shiftpressed',0);
setappdata(resultfig,'cmdpressed',0);
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
if ~strcmp(self.figmode,'lazyupdate')
    tempfile = 'lead_temp_screendump.png';
    self.source_fig = gcf;
    current_paperpositionmode = get(self.source_fig,'PaperPositionMode');
    current_inverthardcopy = get(self.source_fig,'InvertHardcopy');
    set(self.source_fig,'PaperPositionMode','auto');
    set(self.source_fig,'InvertHardcopy','off');
    print(self.source_fig,['-r',num2str(1.5*screen_DPI*self.K(1))], '-dpng', tempfile);
    % change 1.5 to e.g. 2 to get even higher resolution image out.
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
            raw_lowres = double(cat(3,...
                a1(2:self.K(2):end,2:self.K(2):end),...
                a2(2:self.K(2):end,2:self.K(2):end),...
                a3(2:self.K(2):end,2:self.K(2):end)));
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
if strcmp(self.figmode,'figure')
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
elseif strcmp(self.figmode,'publish')
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
elseif strcmp(self.figmode,'update')
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
elseif strcmp(self.figmode,'lazyupdate')
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
