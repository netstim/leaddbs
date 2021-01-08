classdef ea_sweetspot_hemi < handle
    % ROI class to plot niftis on lead dbs resultfig / 3D Matlab figures
    % example:
    % figure; pobj.openedit=1; ea_sweetspot_hemi([spm('dir'),filesep,'toolbox',filesep,'OldSeg',filesep,'grey.nii'],pobj); a=light; axis('off','equal')
    % A. Horn

    properties (SetObservable)
        niftiFilename % original nifti filename
        nii % nifti loaded
        threshold % threshold to visualize
        thresholdType % what the threshold should represent
        color % color of patch
        alpha=0.7 % alpha of patch
        fv % faces and vertices of patch
        sfv % smoothed version
        cdat % color data of patch
        nearelectrodes % handles to ea_trajectory objects plotted in the same figure (by elvis)
        visible='on' % turn on/off
        name % name to be shown
        smooth % smooth by FWHM
        binary % is binary ROI
        hullsimplify % simplify hull
        fmix % factor mix
        peak % peak coordinate
        heatcontacts=1 % heat surrounding electrode contacts
        hail % line to closest contact
        max % maxvalue in nifti - 0.1 of the way
        min % minvalue in nifti + 0.1 of the way
        controlH % handle to color / threshold control figure
        plotFigureH % handle of figure on which to plot
        patchH % handle of patch
        toggleH % toggle handle
        htH % handle for toggle toolbar
    end

    methods
        function obj=ea_sweetspot_hemi(niftiFilename,pobj) % generator function
            if exist('niftiFilename','var') && ~isempty(niftiFilename)
                obj.niftiFilename=niftiFilename;
            end

            try
                obj.name=pobj.name;
            catch
                [~,obj.name]=fileparts(obj.niftiFilename);
            end

            obj.plotFigureH=gcf;

            if exist('pobj','var') && ~isempty(pobj)
                try
                    obj.plotFigureH=pobj.plotFigureH;
                end
            end
            try
                obj.htH=pobj.htH;
            catch
                obj.htH=getappdata(obj.plotFigureH,'addht');
            end
            if isempty(obj.htH) % first ROI
                obj.htH=uitoolbar(obj.plotFigureH);
                setappdata(obj.plotFigureH,'addht',obj.htH);
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

            % load nifti
            obj.nii=load(obj.niftiFilename);
            if ~isfield(obj.nii,'thresholdType')
                obj.nii.thresholdType='absolute';
            end
            if ~isfield(obj.nii,'thresholdFactor')
                obj.nii.thresholdFactor=2;
            end
            obj.binary=0;
            options.prefs=ea_prefs;
            obj.max=ea_nanmax(obj.nii.img(~(obj.nii.img==0)));
            obj.min=ea_nanmin(obj.nii.img(~(obj.nii.img==0)));
            maxmindiff=obj.max-obj.min;
            obj.max=obj.max-0.1*maxmindiff;
            obj.min=obj.min+0.1*maxmindiff;
            obj.fmix=ones(length(obj.nii.factors),1)./length(obj.nii.factors);
            try
                obj.threshold=pobj.threshold;
            catch
                switch obj.nii.thresholdType
                    case 'absolute'
                        if obj.binary
                            obj.threshold=obj.max/2;
                        else
                            obj.threshold=obj.max-0.5*maxmindiff;
                        end
                    case 'percent'
                        obj.threshold=0.95;
                end
            end

            obj.smooth=options.prefs.hullsmooth;
            obj.hullsimplify=options.prefs.hullsimplify;
            set(0,'CurrentFigure',obj.plotFigureH);
            obj.patchH=patch;

            obj.toggleH=uitoggletool;
            % Get the underlying java object using findobj
            jtoggle = findjobj(obj.toggleH);
            % Specify a callback to be triggered on any mouse release event
            set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})

            update_sweetspot(obj);

            addlistener(obj, 'visible', 'PostSet', @obj.changeevent);
            addlistener(obj, 'color', 'PostSet', @obj.changeevent);
            addlistener(obj, 'threshold', 'PostSet', @obj.changeevent);
            addlistener(obj, 'smooth', 'PostSet', @obj.changeevent);
            addlistener(obj, 'hullsimplify', 'PostSet', @obj.changeevent);
            addlistener(obj, 'alpha', 'PostSet', @obj.changeevent);
            addlistener(obj, 'fmix', 'PostSet', @obj.changeevent);

            if exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit
            	ea_editroi([],[],obj)
            end
        end

        function changeevent(~, ~, event)
            update_sweetspot(event.AffectedObject,event.Source.Name);
        end

        function obj=update_sweetspot(obj,evtnm) % update ROI
            if ~exist('evtnm','var')
                evtnm='all';
            end

            obj.nearelectrodes=getappdata(obj.plotFigureH,'el_render');

            if ismember(evtnm,{'all','threshold','smooth','hullsimplify','fmix'}) % need to recalc fv here:
                obj.nii.img(:)=sum(obj.nii.X.*repmat(obj.fmix,1,size(obj.nii.X,2)));
                [~,idx]=ea_nanmax(obj.nii.img(:));
                [x,y,z]=ind2sub(size(obj.nii.img),idx);
                obj.peak=obj.nii.mat*[x;y;z;1];
                obj.peak=obj.peak(1:3);
                bb=[1,1,1;size(obj.nii.img)];
                bb=ea_vox2mm(bb,obj.nii.mat);
                gv=cell(3,1);
                for dim=1:3
                    gv{dim}=linspace(bb(1,dim),bb(2,dim),size(obj.nii.img,dim));
                end
                [X,Y,Z]=meshgrid(gv{1},gv{2},gv{3});

                switch obj.nii.thresholdType
                    case 'absolute' % threshold at absolute percentage
                        obj.fv=isosurface(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.threshold);
                    case 'percent' % threshold by showing percent of peak voxels
                        ints=sort(obj.nii.img(:),'descend');
                        tthresh=ints(round((((obj.threshold)/100)/obj.nii.thresholdFactor)*length(ints)));
                        obj.fv=isosurface(X,Y,Z,permute(obj.nii.img,[2,1,3]),tthresh);
                end
                fvc=isocaps(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.threshold);
                obj.fv.faces=[obj.fv.faces;fvc.faces+size(obj.fv.vertices,1)];
                obj.fv.vertices=[obj.fv.vertices;fvc.vertices];

                if obj.smooth
                    obj.sfv=ea_smoothpatch(obj.fv,1,obj.smooth);
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
                if obj.binary
                    obj.cdat=abs(repmat(atlasc,length(obj.sfv.vertices),1) ... % C-Data for surface
                        +randn(length(obj.sfv.vertices),1)*2)';
                else
                    obj.cdat=isocolors(X,Y,Z,permute(obj.nii.img,[2,1,3]),obj.sfv.vertices);
                end
            end

            jetlist=jet;

            co=ones(1,1,3);
            co(1,1,:)=obj.color;
            atlasc=double(rgb2ind(co,jetlist));

            if obj.heatcontacts
                colorconts(obj)
            end

            % surface patch
            set(0,'CurrentFigure',obj.plotFigureH);
            set(obj.patchH,...
                {'Faces','Vertices','CData','FaceColor','FaceAlpha','EdgeColor','FaceLighting','Visible'},...
                {obj.sfv.faces,obj.sfv.vertices,obj.cdat,obj.color,obj.alpha,'none','phong',obj.visible});

            % toggle button
            set(obj.toggleH,...
                {'Parent','CData','TooltipString','OnCallback','OffCallback','State'},...
                {obj.htH,ea_get_icn('atlas',obj.color),stripext(obj.niftiFilename),{@ea_sweetspotvisible,'on',obj},{@ea_sweetspotvisible,'off',obj},obj.visible});
        end
    end
end


function colorconts(obj)
% now color surrounding electrodes by closeness
    if ~isempty(obj.nearelectrodes)
        els=obj.nearelectrodes;
        % hold on
        % plot3(obj.peak(1),obj.peak(2),obj.peak(3),'y*');
        elrender=getappdata(obj.plotFigureH,'el_render');
        for el=1:length(els)
            cpatches=elrender(el).elpatch(logical(els(el).eltype)); % patches of contacts
            for cont=1:length(els(el).elstruct.coords_mm{els(el).side})
                [~,D(cont)]=knnsearch(cpatches(cont).Vertices,mean(obj.fv.vertices));
            end
            ccol=repmat(0.3,length(els(el).elstruct.coords_mm{els(el).side}),3);

            ofix=find(D<5);
            [~,ix]=sort(D(D<5));

            if ~isempty(ix)
                SweetspotToggleTag = ['SweetspotToggleEle', num2str(el)];
                if isempty(get(obj.toggleH.Parent.Children(1), 'Tag'))
                    set(obj.toggleH.Parent.Children(1), 'Tag', SweetspotToggleTag);
                end

                SweetspotPatchTag = ['SweetspotPatchEle', num2str(el)];
                if isempty(get(obj.patchH.Parent.Children(1), 'Tag'))
                    set(obj.patchH.Parent.Children(1), 'Tag', SweetspotPatchTag);
                end

                allToggle = findobj(obj.toggleH.Parent.Children, 'Tag', SweetspotToggleTag, 'State', 'on');
                allToggleColor = cell2mat(arrayfun(@(x) reshape(x.CData(1,1,:),1,3), allToggle, 'Uniform',0));
                if isempty(allToggle)   % ALL sweetspot toggletools are in 'OFF' state
                    if isempty(obj.toggleH.CData)  % initializing the 1st toggletools w.r.t. one electrode
                        ccol = rgb2ccol(obj.color, ccol, D, ofix, ix);
                    else  % turning OFF all sweetpot toggletools
                        ccol = [];
                    end
                else    % NOT ALL sweetspot toggletools are in 'OFF' state
                    if isempty(obj.toggleH.CData)  % initializing other sweetspots w.r.t. one electrode
                        newColor = mean([allToggleColor; obj.color]);
                        ccol = rgb2ccol(newColor, ccol, D, ofix, ix);
                    elseif strcmp(get(obj.toggleH, 'State'), 'off')  % turning OFF sweetpot toggletool
                        newColor = mean(allToggleColor, 1);
                        ccol = rgb2ccol(newColor, ccol, D, ofix, ix);
                    else  % turning ON sweetpot toggletool
                        newColor = mean(allToggleColor, 1);
                        ccol = rgb2ccol(newColor, ccol, D, ofix, ix);
                    end
                end

                els(el).colorMacroContacts=ccol;

                % draw line to closest contact
                delete(obj.hail);
                if isempty(obj.toggleH.CData) || ...  % initializing sweetspot
                   strcmp(get(obj.toggleH, 'State'), 'on') % turning on sweetspot
                    Pline=[mean(obj.fv.vertices);els(el).elstruct.coords_mm{els(el).side}(ofix(ix(1)),:)];
                    hold on
                    [obj.hail,fv]=ea_plot3t(Pline(:,1),Pline(:,2),Pline(:,3),0.05,[0.1,0.1,0.1],12,1);
                end
            end
        end
    end
end


function rightcallback(src, evnt, obj)
    if evnt.getButton() == 3
        ea_editroi(src,evnt,obj)
    end
end


function ea_editroi(Hobj,evt,obj)
    obj.controlH = ea_sweetspotcontrol_hemi(obj);
end


function ea_sweetspotvisible(Hobj,evt,onoff,obj)
    obj.visible = onoff;
end


function fn = stripext(fn)
    [~,fn] = fileparts(fn);
end


function ccol = rgb2ccol(rgb, ccol, D, ofix, ix)
    P = repmat(rgb,length(ix),1);
    cweights=(1./exp(D(ix))).^2;
    cweights=(cweights/ea_nanmax(cweights))';
    ccol(ofix(ix),:)=ccol(ofix(ix),:).*repmat(1-cweights,1,3)+P.*repmat(cweights,1,3);
end
