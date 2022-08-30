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
    end

    methods

        function sobj = saveobj(obj)
            fn=fieldnames(obj);
            for f=1:length(fn)
                if ~ismember(fn{f},{'plotFigureH','patchH','toggleH','htH','controlH'})
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

                try
                    obj.plotFigureH=pobj.plotFigureH;
                catch
                    obj.plotFigureH=gcf;
                end

                try
                    obj.htH=pobj.htH;
                    if isempty(obj.htH)
                        obj.htH=getappdata(obj.plotFigureH,'addht');
                    end
                catch
                    obj.htH=getappdata(obj.plotFigureH,'addht');
                end

                if isempty(obj.htH) % first ROI
                    obj.htH=uitoolbar(obj.plotFigureH);
                end

                set(0,'CurrentFigure',obj.plotFigureH);
                % set cdata
                if exist('pobj','var') && ~isempty(pobj)
                    try
                        obj.color=pobj.color;
                    catch
                        obj.color = ea_uisetcolor;
                    end
                else
                    obj.color = ea_uisetcolor;
                end

                try
                    obj.binary=pobj.binary;
                end

                try
                    obj.usesolidcolor=pobj.usesolidcolor;
                end

                try
                    obj.colormap=pobj.colormap;
                end

                try
                    obj.nii=pobj.nii;
                catch
                    % load nifti
                    obj.nii=ea_load_nii(obj.niftiFilename);
                    obj.nii.img(obj.nii.img==0) = nan;
                    obj.nii.img(isinf(obj.nii.img)) = nan;

                    if length(unique(obj.nii.img(~isnan(obj.nii.img))))==1
                        obj.binary=1;
                    else
                        obj.nii.img=obj.nii.img-ea_nanmin(obj.nii.img(:)); % set min to zero
                        obj.binary=0;
                    end
                    obj.nii.img(isnan(obj.nii.img)) = 0;
                    obj.nii.img(isinf(obj.nii.img)) = 0;
                end
                options.prefs=ea_prefs;

                obj.max=ea_nanmax(obj.nii.img(~(obj.nii.img==0)));
                obj.min=ea_nanmin(obj.nii.img(~(obj.nii.img==0)));
                maxmindiff=obj.max-obj.min;
                obj.max=obj.max; %-0.1*maxmindiff;
                obj.min=obj.min; %+0.1*maxmindiff;
                try
                    obj.threshold=pobj.threshold;
                catch
                    if obj.binary
                        obj.threshold=obj.max/2;
                    else
                        obj.threshold=obj.max-0.5*maxmindiff;
                    end
                end
                try
                    obj.smooth=pobj.smooth;
                catch
                    obj.smooth=options.prefs.hullsmooth;
                end
                try
                    obj.hullsimplify=pobj.hullsimplify;
                catch
                    obj.hullsimplify=options.prefs.hullsimplify;
                end
                set(0,'CurrentFigure',obj.plotFigureH);
                obj.patchH=patch;

                obj.toggleH=uitoggletool;

                update_roi(obj);
                breathelife(obj);
            end

            if exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit
                ea_editroi([],[],obj)
            end
        end

        function breathelife(obj)
            addlistener(obj,'Visible','PostSet',...
                @changeevent);
            addlistener(obj,'color','PostSet',...
                @changeevent);
            addlistener(obj,'usesolidcolor','PostSet',...
                @changeevent);
            addlistener(obj,'threshold','PostSet',...
                @changeevent);
            addlistener(obj,'smooth','PostSet',...
                @changeevent);
            addlistener(obj,'hullsimplify','PostSet',...
                @changeevent);
            addlistener(obj,'alpha','PostSet',...
                @changeevent);
            addlistener(obj,'SpecularColorReflectance','PostSet',...
                @changeevent);
            addlistener(obj,'SpecularExponent','PostSet',...
                @changeevent);
            addlistener(obj,'SpecularStrength','PostSet',...
                @changeevent);
            addlistener(obj,'DiffuseStrength','PostSet',...
                @changeevent);
            addlistener(obj,'AmbientStrength','PostSet',...
                @changeevent);
            addlistener(obj,'edgecolor','PostSet',...
                @changeevent);

            if isempty(obj.toggleH)
                obj.toggleH=uitoggletool(obj.htH);
            end
            if strcmp(obj.plotFigureH.Visible,'on') % only make GUI functioning if figure is visible.
                % Get the underlying java object using findobj
                jtoggle = findjobj(obj.toggleH);
                % Specify a callback to be triggered on any mouse release event
                set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})
            end
            setappdata(obj.plotFigureH,'addht',obj.htH);
        end

        function obj=update_roi(obj,evtnm) % update ROI
            if ~exist('evtnm','var')
                evtnm='all';
            end
            if isempty(obj.patchH)
                obj.patchH=patch;
            end
            if ismember(evtnm,{'all','threshold','smooth','hullsimplify','usesolidcolor'}) % need to recalc fv here:
                bb = [1,1,1;size(obj.nii.img)];
                bb = ea_vox2mm(bb, obj.nii.mat);
                gv=cell(3,1);
                for dim=1:3
                    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(obj.nii.img,dim));
                end
                [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});

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

                if obj.binary || obj.usesolidcolor
                    obj.cdat=abs(repmat(obj.color,length(obj.sfv.vertices),1) ... % C-Data for surface
                        +randn(length(obj.sfv.vertices),1)*2)';
                else
                    if isempty(obj.colormap) % make sure some colormap is set.
                        obj.colormap = ea_colorgradient(length(gray), [0,0,1], [1,1,1], [1,0,0]); % default blue to red colormap
                    end
                    obj.cdat=isocolors(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.sfv.vertices);
                    obj.cdat=round((ea_contrast(obj.cdat).*(length(gray)-1))+1);
                    obj.cdat=obj.colormap(obj.cdat,:);
                end
            end

            %co=ones(1,1,3);
            %co(1,1,:)=obj.color;
            %atlasc=double(rgb2ind(co,jetlist));

            % show atlas.
            set(0,'CurrentFigure',obj.plotFigureH);

            if isempty(obj.edgecolor)
                obj.edgecolor='none';
            end

            % Set obj Tag
            if ~isempty(obj.Tag)
                roiTag = obj.Tag;
            else
                roiTag = obj.name;
            end

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

            % add toggle button:
            set(obj.toggleH,...
                {'Parent','CData','TooltipString','OnCallback','OffCallback','State','Tooltip'},...
                {obj.htH,ea_get_icn('atlas',obj.color),stripext(obj.niftiFilename),{@ea_roivisible,'on',obj},{@ea_roivisible,'off',obj},obj.Visible,roiTag});
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

function ea_editroi(Hobj,evt,obj)
    obj.controlH=ea_roicontrol(obj);
end

function ea_roivisible(Hobj,evt,onoff,obj)
    obj.Visible=onoff;
end

function fn=stripext(fn)
    [~,fn]=fileparts(fn);
end
