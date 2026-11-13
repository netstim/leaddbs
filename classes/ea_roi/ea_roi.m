classdef ea_roi < handle
    % ROI class to plot niftis on lead dbs resultfig / 3D Matlab figures
    % example:
    % figure; pobj.openedit=1; ea_roi([spm('dir'),filesep,'toolbox',filesep,'OldSeg',filesep,'grey.nii'],pobj); a=light; axis('off','equal')
    % A. Horn

    properties (SetObservable)
        niftiFilename % original nifti filename
        nii % nifti loaded
        threshold % threshold to visualize
        color % facecolor of patch
        edgecolor='none' % edgecolor of patch
        usesolidcolor=1 % whether to use isocolors or a solid manually defined color
        alpha=0.7 % alpha of patch
        fv % faces and vertices of patch
        sfv % smoothed version
        cdat % color data of patch
        Visible='on' % turn on/off
        name % name to be shown
        smooth % smooth by FWHM
        binary % is binary ROI
        hullsimplify % simplify hull
        max % maxvalue in nifti
        min % minvalue in nifti
        controlH % handle to color / threshold control figure
        plotFigureH % handle of figure on which to plot
        patchH % handle of patch
        colormap % for nonbinary ROI
        toggleH % toggle handle
        htH % handle for toggle toolbar
        Tag % tag of ROI can be used in multi-roi scenes
        SpecularColorReflectance=1 % patch property
        SpecularExponent=3 % patch property
        SpecularStrength=0.3 % patch property
        DiffuseStrength=0.4 % patch property
        AmbientStrength=0.3 % patch property
    end

    properties (Access = private, Transient)
        listenersInitialized = false
        contextToggleH = []
    end

    methods(Static)

        function obj=loadobj(s)
            if isstruct(s)
                newObj=ea_roi();
                fn=fieldnames(s);
                for f=1:length(fn)
                    newObj.(fn{f})=s.(fn{f});
                end
                obj=newObj;
            else
                obj=s;
            end
        end

        function bindcontextmenus(rois,resultfig)
            % Bind atlas toolbar right clicks in one Java lookup. findjobj
            % flushes the renderer on every call, which is costly when a
            % detailed atlas is still being assembled.
            if ~isMatlabVer('<',[25,1]) || ~strcmp(resultfig.Visible,'on')
                return
            end
            candidates = {};
            for r = 1:numel(rois)
                obj = rois{r};
                if isa(obj,'ea_roi') && ~isempty(obj.toggleH) && ...
                        ~isequal(obj.contextToggleH,obj.toggleH)
                    candidates{end+1} = obj; %#ok<AGROW>
                end
            end
            if isempty(candidates)
                return
            end

            tooltips = cell(size(candidates));
            for r = 1:numel(candidates)
                tooltips{r} = get(candidates{r}.toggleH,'TooltipString');
            end
            cleanup = onCleanup(@() restoretooltips(candidates,tooltips));
            for r = 1:numel(candidates)
                % Temporary unique tooltips also distinguish repeated atlas
                % names and left/right structures with identical labels.
                set(candidates{r}.toggleH,'TooltipString',sprintf('ea_roi_context_%d',r));
            end
            try
                jtoggles = findjobj(resultfig,'class','ToggleButton','nomenu');
                for j = 1:numel(jtoggles)
                    tooltip = char(jtoggles(j).getToolTipText);
                    r = sscanf(tooltip,'ea_roi_context_%d');
                    if isscalar(r) && r >= 1 && r <= numel(candidates) && ...
                            strcmp(tooltip,sprintf('ea_roi_context_%d',r))
                        obj = candidates{r};
                        set(jtoggles(j),'MouseReleasedCallback',{@rightcallback,obj});
                        jtoggles(j).setToolTipText(tooltips{r});
                        obj.contextToggleH = obj.toggleH;
                    end
                end
            catch
                % Retain the established per-button lookup as a fallback
                % for Java toolbar implementations with a different tree.
            end
            clear cleanup
            for r = 1:numel(candidates)
                if ~isequal(candidates{r}.contextToggleH,candidates{r}.toggleH)
                    candidates{r}.breathelife;
                end
            end
        end
    end

    methods

        function sobj = saveobj(obj)
            fn=fieldnames(obj);
            for f=1:length(fn)
                if ~ismember(fn{f},{'plotFigureH','patchH','toggleH','htH','controlH',...
                        'listenersInitialized','contextToggleH'})
                    sobj.(fn{f})=obj.(fn{f});
                end
            end
        end

        function obj=ea_roi(niftiFilename,pobj) % generator function
            if nargin
                if exist('niftiFilename','var') && ~isempty(niftiFilename)
                    obj.niftiFilename=niftiFilename;
                end

                try
                    obj.name=pobj.name;
                catch
                    [~,obj.name]=ea_niifileparts(obj.niftiFilename);
                end
                
                obj.Tag = obj.name;

                if isfield(pobj, 'plotFigureH')
                    obj.plotFigureH = pobj.plotFigureH;
                else
                    currentfig = gcf;
                    stack = dbstack;
                    if contains(currentfig.Name, 'Electrode-Scene') && ~ismember('ea_genatlastable', {stack.name})
                        % Set figure and toolbar only in elvis
                        obj.plotFigureH = currentfig;
                    end
                end

                if isfield(pobj, 'htH')
                    obj.htH = pobj.htH;
                end

                if isempty(obj.htH) && ~isempty(obj.plotFigureH)
                    obj.htH = getappdata(obj.plotFigureH, 'addht');
                    if isempty(obj.htH)
                        obj.htH = uitoolbar(obj.plotFigureH);  % first ROI
                    end
                end

                if ~isempty(obj.plotFigureH)
                    set(0,'CurrentFigure',obj.plotFigureH);
                end

                % set cdata
                if exist('pobj','var') && ~isempty(pobj)
                    if isfield(pobj, 'color')
                        obj.color = pobj.color;
                    else
                        obj.color = ea_uisetcolor;
                    end
                else
                    obj.color = ea_uisetcolor;
                end

                if isfield(pobj, 'binary')
                    obj.binary = pobj.binary;
                end

                if isfield(pobj, 'usesolidcolor')
                    obj.usesolidcolor = pobj.usesolidcolor;
                end

                if isfield(pobj, 'colormap')
                    obj.colormap = pobj.colormap;
                end

                if isfield(pobj, 'nii')
                    obj.nii = pobj.nii;
                else
                    % load nifti
                    obj.nii=ea_load_nii(obj.niftiFilename);
                end

                obj.nii.img(obj.nii.img==0) = nan;
                obj.nii.img(isinf(obj.nii.img)) = nan;

                if isscalar(unique(obj.nii.img(~isnan(obj.nii.img))))
                    obj.binary=1;
                else
                    obj.nii.img=obj.nii.img-ea_nanmin(obj.nii.img(:)); % set min to zero
                    obj.binary=0;
                end
                obj.nii.img(isnan(obj.nii.img)) = 0;
                obj.nii.img(isinf(obj.nii.img)) = 0;
                options.prefs=ea_prefs;

                obj.max=ea_nanmax(obj.nii.img(~(obj.nii.img==0)));
                obj.min=ea_nanmin(obj.nii.img(~(obj.nii.img==0)));
                maxmindiff=obj.max-obj.min;
                obj.max=obj.max; %-0.1*maxmindiff;
                obj.min=obj.min; %+0.1*maxmindiff;

                if isfield(pobj, 'threshold')
                    obj.threshold = pobj.threshold;
                else
                    if obj.binary
                        obj.threshold = obj.max/2;
                    else
                        obj.threshold = obj.max-0.5*maxmindiff;
                    end
                end

                if isfield(pobj, 'smooth')
                    obj.smooth = pobj.smooth;
                else
                    obj.smooth = options.prefs.hullsmooth;
                end

                if isfield(pobj, 'hullsimplify')
                    obj.hullsimplify = pobj.hullsimplify;
                else
                    obj.hullsimplify = options.prefs.hullsimplify;
                end

                if ~isempty(obj.plotFigureH)
                    set(0,'CurrentFigure',obj.plotFigureH);
                    obj.patchH = patch;
                end

                if ~isempty(obj.htH)
                    obj.toggleH = uitoggletool(obj.htH);
                    stack = dbstack;
                    if ismember('ea_addobj', {stack.name})
                        obj.toggleH.UserData = 'roi';
                    end
                end

                update_roi(obj);
                breathelife(obj);
            end

            if exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit
                ea_editroi([],[],obj)
            end
        end

        function breathelife(obj,bindJava)
            if nargin < 2
                bindJava = true;
            end
            % Loaded atlas ROIs may be attached more than once. A second
            % listener would repeat every geometry rebuild and UI update.
            if ~obj.listenersInitialized
                observed = {'Visible','color','usesolidcolor','threshold',...
                    'smooth','hullsimplify','alpha','SpecularColorReflectance',...
                    'SpecularExponent','SpecularStrength','DiffuseStrength',...
                    'AmbientStrength','edgecolor'};
                for p = 1:numel(observed)
                    addlistener(obj,observed{p},'PostSet',@changeevent);
                end
                obj.listenersInitialized = true;
            end

            if isempty(obj.toggleH) && ~isempty(obj.htH)
                obj.toggleH = uitoggletool(obj.htH);
            end

            if ~isempty(obj.plotFigureH)
                if strcmp(obj.plotFigureH.Visible,'on') % only make GUI functioning if figure is visible.
                    if bindJava && isMatlabVer('<', [25,1]) && ~isequal(obj.contextToggleH,obj.toggleH)
                        % Get the underlying java object using findobj
                        jtoggle = findjobj(obj.toggleH);
                        % Specify a callback to be triggered on any mouse release event
                        set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})
                        obj.contextToggleH = obj.toggleH;
                    end
                end
                setappdata(obj.plotFigureH, 'addht', obj.htH);
            end
        end

        function [obj,geometryRebuilt]=update_roi(obj,evtnm) % update ROI
            geometryRebuilt = false;
            if ~exist('evtnm','var')
                evtnm='all';
            end

            if isempty(obj.patchH) && ~isempty(obj.plotFigureH)
                obj.patchH = patch;
            end
            if isempty(obj.toggleH) && ~isempty(obj.htH)
                obj.toggleH = uitoggletool(obj.htH);
            end

            % Visibility, appearance and toolbar state do not change the
            % mesh. In particular, avoid resending large Faces/Vertices
            % arrays to the renderer when checking an atlas item on/off.
            if strcmp(evtnm,'Visible')
                if ~isempty(obj.patchH)
                    set(obj.patchH,'Visible',obj.Visible);
                end
                if ~isempty(obj.toggleH)
                    set(obj.toggleH,'State',obj.Visible);
                end
                return
            elseif ismember(evtnm,{'alpha','edgecolor','SpecularColorReflectance',...
                    'SpecularExponent','SpecularStrength','DiffuseStrength','AmbientStrength'})
                if ~isempty(obj.patchH)
                    patchProperty = evtnm;
                    if strcmp(evtnm,'alpha')
                        patchProperty = 'FaceAlpha';
                    elseif strcmp(evtnm,'edgecolor') && isempty(obj.edgecolor)
                        obj.edgecolor = 'none';
                    end
                    set(obj.patchH,patchProperty,obj.(evtnm));
                end
                return
            elseif strcmp(evtnm,'color')
                if ~isempty(obj.patchH) && (obj.binary || obj.usesolidcolor)
                    set(obj.patchH,'FaceColor',obj.color);
                end
                if ~isempty(obj.toggleH)
                    set(obj.toggleH,'CData',ea_get_icn('atlas',obj.color));
                end
                return
            end

            rebuildGeometry = ismember(evtnm,{'all','threshold','smooth','hullsimplify'});
            if strcmp(evtnm,'initialize')
                % Atlas indexes already contain both original and displayed
                % meshes. The caller requests 'all' if mesh settings changed.
                rebuildGeometry = ~isstruct(obj.fv) || ~isstruct(obj.sfv) || ...
                    ~isfield(obj.fv,'vertices') || ~isfield(obj.fv,'faces') || ...
                    ~isfield(obj.sfv,'vertices') || ~isfield(obj.sfv,'faces') || ...
                    isempty(obj.sfv.vertices) || isempty(obj.sfv.faces);
            end
            geometryRebuilt = rebuildGeometry;
            refreshColors = rebuildGeometry || strcmp(evtnm,'usesolidcolor') || ...
                (strcmp(evtnm,'initialize') && ~(obj.binary || obj.usesolidcolor) && ...
                size(obj.cdat,1) ~= size(obj.sfv.vertices,1));

            if rebuildGeometry || (refreshColors && ~(obj.binary || obj.usesolidcolor))
                img_size = size(obj.nii.img);
                [X_vox, Y_vox, Z_vox] = meshgrid(1:img_size(1), 1:img_size(2), 1:img_size(3));
                XYZ_mm = ea_vox2mm([X_vox(:), Y_vox(:), Z_vox(:)], obj.nii.mat);
                X = reshape(XYZ_mm(:,1), size(X_vox));
                Y = reshape(XYZ_mm(:,2), size(Y_vox));
                Z = reshape(XYZ_mm(:,3), size(Z_vox));
            end

            if rebuildGeometry
                obj.fv=isosurface(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.threshold);
                fvc=isocaps(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.threshold);
                obj.fv.faces=[obj.fv.faces;fvc.faces+size(obj.fv.vertices,1)];
                obj.fv.vertices=[obj.fv.vertices;fvc.vertices];

                if obj.smooth
                    if ~isempty(obj.fv.vertices) && ~isempty(obj.fv.faces)
                        
                        obj.sfv=ea_smoothpatch(obj.fv,1,obj.smooth);
                    else
                        return
                    end
                else
                    obj.sfv=obj.fv;
                end

                if ischar(obj.hullsimplify)
                    % get to 700 faces
                    simplify=700/length(obj.sfv.faces);
                    obj.sfv=reducepatch(obj.sfv,simplify);
                else
                    if obj.hullsimplify<1 && obj.hullsimplify>0
                        obj.sfv=reducepatch(obj.sfv,obj.hullsimplify);
                    elseif obj.hullsimplify>1
                        simplify=obj.hullsimplify/length(obj.fv.faces);
                        obj.sfv=reducepatch(obj.sfv,simplify);
                    end
                end

            end

            if refreshColors
                if obj.binary || obj.usesolidcolor
                    obj.cdat=abs(repmat(obj.color,length(obj.sfv.vertices),1) ... % C-Data for surface
                        +randn(length(obj.sfv.vertices),1)*2)';
                else
                    if isempty(obj.colormap) % make sure some colormap is set.
                        obj.colormap = ea_colorgradient(length(gray), [0,0,1], [1,1,1], [1,0,0]); % default blue to red colormap
                    end
                    obj.cdat=isocolors(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.sfv.vertices);
                    ncolors=size(obj.colormap,1);
                    obj.cdat=round((ea_contrast(obj.cdat).*(ncolors-1))+1);
                    obj.cdat(~isfinite(obj.cdat))=1;
                    obj.cdat=max(1,min(ncolors,obj.cdat));
                    obj.cdat=obj.colormap(obj.cdat,:);
                end
            end

            if strcmp(evtnm,'usesolidcolor')
                if ~isempty(obj.patchH)
                    if obj.binary || obj.usesolidcolor
                        set(obj.patchH,'FaceColor',obj.color);
                    else
                        set(obj.patchH,'FaceVertexCData',obj.cdat,'FaceColor','interp');
                    end
                end
                return
            end

            %co=ones(1,1,3);
            %co(1,1,:)=obj.color;
            %atlasc=double(rgb2ind(co,jetlist));

            if isempty(obj.edgecolor)
                obj.edgecolor='none';
            end

            % Set obj Tag
            if ~isempty(obj.Tag)
                roiTag = obj.Tag;
            else
                roiTag = obj.name;
            end

            % show atlas.
            if ~isempty(obj.plotFigureH)
                set(0,'CurrentFigure',obj.plotFigureH);
                set(obj.patchH,...
                    {'Faces','Vertices','FaceAlpha','EdgeColor','EdgeLighting','FaceLighting','Visible','SpecularColorReflectance','SpecularExponent','SpecularStrength','DiffuseStrength','AmbientStrength','Tag'},...
                    {obj.sfv.faces,obj.sfv.vertices,obj.alpha,obj.edgecolor,'gouraud','gouraud',obj.Visible,obj.SpecularColorReflectance,obj.SpecularExponent,obj.SpecularStrength,obj.DiffuseStrength,obj.AmbientStrength,roiTag});
                if obj.binary || obj.usesolidcolor
                    set(obj.patchH,...
                        {'FaceColor'},...
                        {obj.color});
                else
                    set(obj.patchH,...
                        {'FaceVertexCData','FaceColor'},...
                        {obj.cdat,'interp'});
                end
            end

            % add toggle button:
            if ~isempty(obj.htH)
                set(obj.toggleH,...
                    {'Parent','CData','TooltipString','OnCallback','OffCallback','State','Tag','Tooltip'},...
                    {obj.htH,ea_get_icn('atlas',obj.color),stripext(obj.niftiFilename),{@ea_roivisible,'on',obj},{@ea_roivisible,'off',obj},obj.Visible,roiTag,roiTag});
            end
        end

        function delete(obj)
            delete(obj.patchH);
        end
    end
end

function changeevent(~,event)
    update_roi(event.AffectedObject,event.Source.Name);
end

function rightcallback(src, evnt,obj)
    if evnt.getButton() == 3
        ea_editroi(src,evnt,obj)
    end
end

function ea_editroi(~,~,obj)
    obj.controlH = ea_roicontrol(obj);
end

function ea_roivisible(~,~,onoff,obj)
    obj.Visible=onoff;
end

function restoretooltips(rois,tooltips)
for r = 1:numel(rois)
    if isgraphics(rois{r}.toggleH)
        set(rois{r}.toggleH,'TooltipString',tooltips{r});
    end
end
end

function fn=stripext(fn)
    [~,fn]=fileparts(fn);
end
