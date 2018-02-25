classdef ea_trajectory < handle
    % Trajectory class to handle trajectories in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        elstruct % reconstruction of electrodes as handled by ea_elvis
        elpatch % handle to macroelectrode patch
        eltype % indexes 1 for electrode contacts in elpatch
        ellabel % handle to electrode label
        elmodel % elmodel to display
        side % right hemisphere=1, left=2, further sites planned to be possible in the future
        target % target and entrypoints as used in surgical planning
        alpha=0.7 % alpha of Planning patch
        radius=0.2 % radius of Planning line
        color=[0.8,0.3,0.2] % color of Planning patch
        colorMacroContacts=[] % optional coloring of macroelectrode contacts
        options % lead-dbs options struct
        planRelative=[2,1,1,1,1] % First entry: AC=1, MCP=2, PC=3; Second entry: Right=1, Left=2; Third entry: Anterior=1, Posterior=2; Fourth entry: Ventral=1; Dorsal=2; Last entry: ACPC=1, native=2, MNI/Template=3
        hasPlanning % determines if object has information to show a fiducial
        hasMacro % determines if object has information to show a macroelectrode
        relateMicro='macro' % determines if microelectrodes shown should be related to planning Fiducial ('planning') or Macroelectrodes ('macro')
        showPlanning=1 % show planning fiducial
        showMacro=0 % show definitive DBS / macro electrode
        showMicro=0 % show microelectrodes
        merstruct % MERstate object of mer-fiducials
        controlH % handle to trajectory control figure
        plotFigureH % handle of figure on which to plot
        patchMacro % handle of macroelectrode patch
        patchPlanning % handle of planning fiducial patch
        patchMicro % handle of microelectrodes
        toggleH % togglebutton handle that will open planning fiducial control
        htH % handle for toggle toolbar on which toggleH is displayed
        togglestates % show/hide states of primitive toggle button
        toggledefault % which part to show by activating toggletool if none is shown
    end

    methods
        function obj=ea_trajectory(pobj) % class constructor
            try
                obj.plotFigureH=pobj.plotFigureH;
            catch
                obj.plotFigureH=gcf;
            end

            obj.htH=getappdata(obj.plotFigureH,'ht');

            if isempty(obj.htH) % first Entry on toolbar
                obj.htH=uitoolbar(obj.plotFigureH);
                setappdata(obj.plotFigureH,'addht',obj.htH);
            end

            set(0,'CurrentFigure',obj.plotFigureH);
            % set cdata
            try
                obj.color=pobj.color;
            end

            %% initialize reco and controlling entries
            try % target
                obj.elstruct=pobj.elstruct;
            catch
                obj.elstruct=struct;
            end

            if ~exist('pobj','var') % create blank trajectory with planning fiducial only
                obj.hasPlanning=1;
                obj.hasMacro=0;
            else % determine if fiducial and macro information is available
                obj.hasMacro=~isempty(obj.elstruct);
                obj.hasPlanning=~isempty(obj.target);
                obj.showMacro=obj.hasMacro;
                obj.showPlanning=obj.hasPlanning*(~obj.showMacro);
            end

            try
                obj.side=pobj.side;
            catch
                obj.side=1;
            end

            %% initialize further content fields based on given struct if given or else empty / random vars
            if obj.hasPlanning
                try % target
                    obj.target=pobj.target;
                catch
                    obj.target=ea_getstandardtarget(1);
                end
            else
                obj.target=struct;
            end

            try % patchFiducial
                obj.patchPlanning=pobj.patchPlanning;
            catch
                obj.patchPlanning=patch('Visible','off');
            end

            try % patchMacro
                obj.patchMacro=pobj.patchMacro;
            catch
                obj.patchMacro=patch('Visible','off');
            end

            try % patchMicro
                obj.patchMicro=pobj.patchMicro;
            catch
                obj.patchMicro.cent=patch('Visible','off');
                obj.patchMicro.ant=patch('Visible','off');
                obj.patchMicro.post=patch('Visible','off');
                obj.patchMicro.lat=patch('Visible','off');
                obj.patchMicro.med=patch('Visible','off');
                obj.patchMicro.marks=struct;
            end

            if obj.hasMacro && ~obj.hasPlanning
                obj.relateMicro='macro'; % relate microelectrodes to macroelectrode
            elseif ~obj.hasMacro && obj.hasPlanning
                obj.relateMicro='planning'; % relate microelectrodes to planning trajectory
            end

            set(0,'CurrentFigure',obj.plotFigureH);

            try
                obj.options=pobj.options;
            catch
                obj.options=getappdata(obj.plotFigureH,'options');
            end

            switch obj.options.leadprod
                case {'dbs','group'}
                    obj.toggledefault='macro';
                case 'or'
                    obj.toggledefault='planning';
            end

            if isempty(obj.elmodel)
                try
                    obj.elmodel=obj.options.elmodel;
                catch
                    elms=ea_resolve_elspec;
                    obj.elmodel=elms{1};
                end
            end

            if isempty(obj.options)
                ea_warning('Patient information not available.');
            end

            obj.toggleH=uitoggletool;
            update_trajectory(obj);

            % Get the underlying java object using findobj
            jtoggle = findjobj(obj.toggleH);

            % Specify a callback to be triggered on any mouse release event
            set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})
            addlistener(obj, 'showPlanning', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'hasPlanning', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'elmodel', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'showMacro', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'showMicro', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'relateMicro', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'color', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'target', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'alpha', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'target', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'planRelative', 'PostSet', @ea_trajectory.changeevent);
            addlistener(obj, 'colorMacroContacts', 'PostSet', @ea_trajectory.changeevent);

            if (exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit) || ~exist('pobj','var')
                ea_trajectorycontrol(obj)
            end
        end
    end

    methods (Static)
        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
    end
end


function obj=update_trajectory(obj,evtnm) % update ROI
    if ~exist('evtnm','var')
        evtnm='all';
    end
    set(0,'CurrentFigure',obj.plotFigureH);
    if ismember(evtnm,{'all','target','reco','planRelative','hasPlanning','showMicro','relateMicro'}) % need to redraw planning fiducials:
        % planning fiducial
        if obj.hasPlanning
            coords=ea_convertfiducials(obj,[obj.target.target;obj.target.entry]);
            tgt=coords(1,:); ent=coords(2,:);
            for dim=1:3
                traj(:,dim)=linspace(ent(dim),tgt(dim),10);
            end
            delete(obj.patchPlanning);
            obj.patchPlanning=ea_plot3t(traj(:,1),traj(:,2),traj(:,3),obj.radius,obj.color,12,1);
        end

        if obj.showMicro
            merviewer=ea_openmerviewer([],[],obj);
            setappdata(obj.plotFigureH,'merviewer',merviewer);
            obj.merstruct=getappdata(merviewer, 'merstruct'); % bind MERstate object to trajectory object.
        else
            try % close potentially existing MER viewer figure.
                close(getappdata(obj.plotFigureH,'merviewer'));
            end
        end
    end

    if ismember(evtnm,{'color'}) % simply change color of patch
        obj.patchPlanning.FaceVertexCData=repmat(obj.color,size(obj.patchPlanning.FaceVertexCData,1),1);
    end

    if ismember(evtnm,{'showPlanning'}) && obj.hasPlanning
        obj.patchPlanning.Visible=ea_bool2onoff(obj.showPlanning);
    end

    if ismember(evtnm,{'all','elmodel','colorMacroContacts'})
        if obj.showMacro
            try
                delete(obj.elpatch{1}{1});
                delete(obj.ellabel(1));
            end
            poptions=obj.options;
            poptions.elmodel=obj.elmodel;
            obj.elstruct.elmodel=obj.elmodel;
            poptions.sides=obj.side;
            poptions.colorMacroContacts=obj.colorMacroContacts;
            el_render=getappdata(obj.plotFigureH,'el_render');

            [obj.elpatch{1},obj.ellabel(1),obj.eltype]=ea_showelectrode(obj.plotFigureH,obj.elstruct,1,poptions);
            if isempty(el_render)
                clear el_render
            end
            el_render(obj.side)=obj.elpatch{1};
            setappdata(obj.plotFigureH,'el_render',el_render);
            if ~isnan(obj.ellabel(1))
                set(obj.ellabel(1),'Visible','off');
            end
        end
    end

    if ismember(evtnm,{'showMacro'}) && obj.hasMacro
        ea_elvisible([],[],obj.elpatch,1,obj.side,ea_bool2onoff(obj.showMacro),obj.options);
    end

    % add toggle button:
    [~,ptname]=fileparts(fileparts(obj.options.root));
    set(obj.toggleH, {'Parent','CData','TooltipString','OnCallback','OffCallback','State'},...
        {obj.htH,ea_get_icn('electrode'), [ptname,' (',num2str(obj.side),')'], ...
        {@ea_trajvisible,'on',obj}, {@ea_trajvisible,'off',obj}, ...
        ea_bool2onoff(any([obj.showPlanning,obj.showMacro,obj.showMicro]))});
end


function ccoords=ea_convertfiducials(obj,coords)
    for coord=1:size(coords,1)
        thiscoord=coords(coord,:);
        switch obj.planRelative(5)
            case 1 % planning in AC/PC system
                cfg.mapmethod=0;
                cfg.acmcpc=obj.planRelative(1);
                cfg.xmm=thiscoord(1); cfg.ymm=thiscoord(2); cfg.zmm=thiscoord(3);
                if obj.planRelative(2)==2; cfg.xmm=-cfg.xmm; end
                if obj.planRelative(3)==2; cfg.ymm=-cfg.ymm; end
                if obj.planRelative(4)==1; cfg.zmm=-cfg.zmm; end

                switch obj.options.native
                    case 1 % need to convert from AC/PC to native
                        cfg.native=1;
                        wp=ea_acpc2mni(cfg,{[obj.options.root,obj.options.patientname,filesep]});
                        ccoords(coord,:)=wp.WarpedPointNative;
                    case 0 % need to convert from AC/PC to template
                        wp=ea_acpc2mni(cfg,{[obj.options.root,obj.options.patientname,filesep]});
                        ccoords(coord,:)=wp.WarpedPointMNI;
                end
            case 2 % planning in native space
                switch obj.options.native
                    case 1 % leave coords as they are
                        ccoords(coord,:)=coords(coord,:);
                    case 0 % need to convert from native to MNI
                        V = ea_open_vol([obj.options.root,obj.options.patientname,filesep,obj.options.prefs.prenii_unnormalized]);
                        thiscoordvox=V.mat\[thiscoord,1]';
                        ccoords(coord,:)=ea_map_coords(thiscoordvox,...
                            [obj.options.root,obj.options.patientname,filesep,obj.options.prefs.prenii_unnormalized],...
                            [obj.options.root,obj.options.patientname,filesep,'y_ea_inv_normparams.nii'], ...
                            [ea_space,obj.options.primarytemplate,'.nii'])';
                end
            case 3 % planning in template space
                switch obj.options.native
                    case 1 % need to convert from MNI to native
                        % from MNI mm to MNI vox:
                        V = ea_open_vol([ea_space,obj.options.primarytemplate,'.nii']);
                        thiscoordvox=V.mat\[thiscoord,1]';
                        ccoords(coord,:)=ea_map_coords(thiscoordvox,...
                            [ea_space,obj.options.primarytemplate,'.nii'],...
                            [obj.options.root,obj.options.patientname,filesep,'y_ea_normparams.nii'], ...
                            [obj.options.root,obj.options.patientname,filesep,obj.options.prefs.prenii_unnormalized])';
                    case 0 % leave coords as they are
                        ccoords(coord,:)=coords(coord,:);
                end
        end
    end
end


function ea_roivisible(Hobj,evt,onoff,obj)
    obj.visible=onoff;
end


function coords=map_coords_proxy(XYZ,V)
    XYZ=[XYZ';ones(1,size(XYZ,1))];

    coords=V.mat*XYZ;
    coords=coords(1:3,:)';
end


function fn=stripext(fn)
    [~,fn]=fileparts(fn);
end


function rightcallback(src, evnt, obj)
    if evnt.getButton() == 3
        ea_editfiducial(src,evnt, obj)
    end
end


function ea_editfiducial(Hobj, evt, obj)
    obj.controlH=ea_trajectorycontrol(obj);
end


function ea_trajvisible(~, ~, onoff, obj)
    if getappdata(obj.plotFigureH,'altpressed') % hide all
        el_render=getappdata(obj.plotFigureH,'el_render');

        for el=1:length(el_render)
            el_render(el).showPlanning=ea_bool2onoff('off');
            el_render(el).showMacro=ea_bool2onoff(onoff);
            el_render(el).showMicro=ea_bool2onoff('off');
        end
    else
        if ishandle(obj.controlH) % control figure exists
            chandles=getappdata(obj.controlH,'chandles');
            obj.showPlanning=get(chandles.showPlanning,'Value');
            obj.showMacro=get(chandles.showMacro,'Value');
        else
            switch onoff
                case 'off'
                    obj.showPlanning=0;
                    obj.showMacro=0;
                    obj.showMicro=0;
                case 'on'
                    obj.showMacro=1;
            end
        end
    end
end
