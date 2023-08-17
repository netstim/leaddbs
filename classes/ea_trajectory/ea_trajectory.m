classdef ea_trajectory < handle
    % Trajectory class to handle trajectories in lead dbs resultfig / 3D Matlab figures
    % A. Horn

    properties (SetObservable)
        elstruct % reconstruction of electrodes as handled by ea_elvis
        plan2elstruct % reconstruction (pseudo) of stereotactical plan
        plan2elstruct_model='Medtronic 3389' % electrode model of pseudo reconstruction of stereotactical plan
        electrodeRelativeToPlan=3;
        elpatch % handle to DBS electrode patch
        eltype % indexes 1 for electrode contacts in elpatch
        ellabel % handle to electrode label
        elmodel % elmodel to display
        side % right hemisphere=1, left=2, further sites planned to be possible in the future
        target % target and entrypoints as used in surgical planning
        alpha=0.7 % alpha of Planning patch
        radius=0.2 % radius of Planning line
        color=[0.5,0.5,0.5] % color of Planning patch
        colorMacroContacts=[] % optional coloring of macroelectrode contacts
        options % lead-dbs options struct
        planRelative=[2,1,1,1,3] % First entry: AC=1, MCP=2, PC=3; Second entry: Right=1, Left=2; Third entry: Anterior=1, Posterior=2; Fourth entry: Ventral=1; Dorsal=2; Last entry: ACPC=1, native=2, MNI/Template=3
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
        pt=1 % used for patient index count in lead group.
        planningAppearance='line' % can be set to 'electrode' to show 3D electrode instead
    end

    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
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
                try
                    obj.hasMacro=pobj.hasMacro;
                catch
                    obj.hasMacro=~isempty(obj.elstruct);
                end
                try
                    obj.hasPlanning=pobj.hasPlanning;
                catch
                    obj.hasPlanning=~isempty(obj.target);
                end
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

            try % patient index when calling from lead group.
                obj.pt=pobj.pt;
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
                    obj.elmodel=obj.elstruct.elmodel;
                catch
                    try
                        obj.elmodel=obj.options.elmodel;
                    catch
                        elms=ea_resolve_elspec;
                        obj.elmodel=elms{1};
                    end
                end
            end

            if isempty(obj.options)
                ea_warning('Patient information not available.');
            end

            obj.toggleH=uitoggletool;
            update_trajectory(obj);

            % Get the underlying java object using findobj
            if strcmp(obj.plotFigureH.Visible,'on') % only needed if figure is really visible (when calling from lead group it could be hidden).
                jtoggle = findjobj(obj.toggleH);

                % Specify a callback to be triggered on any mouse release event
                set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})
                addlistener(obj, 'showPlanning', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'hasPlanning', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'elmodel', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'elstruct', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'showMacro', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'showMicro', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'relateMicro', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'color', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'target', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'alpha', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'target', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'planRelative', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'colorMacroContacts', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'planningAppearance', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'plan2elstruct_model', 'PostSet', @ea_trajectory.changeevent);
                addlistener(obj, 'electrodeRelativeToPlan', 'PostSet', @ea_trajectory.changeevent);
            end
            if (exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit) || ~exist('pobj','var')
                obj.controlH = ea_trajectorycontrol(obj);
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
    ea_busyaction('on',obj.controlH,'trajectory');
    if strcmp(evtnm,'planRelative')
        fromspace=obj.switchedFromSpace;
        tospace=obj.planRelative(5);
        obj.switchedFromSpace=obj.planRelative(5);
        % target
        tcfg.xmm=obj.target.target(1);
        tcfg.ymm=obj.target.target(2);
        tcfg.zmm=obj.target.target(3);
        tcfg.mapmethod=0;
        tcfg.acmcpc=obj.planRelative(1);
        % entry
        ecfg.xmm=obj.target.entry(1);
        ecfg.ymm=obj.target.entry(2);
        ecfg.zmm=obj.target.entry(3);
        ecfg.mapmethod=0;
        ecfg.acmcpc=obj.planRelative(1);
        if fromspace==1 % from ACPC

            if obj.planRelative(2)==2
                tcfg.xmm=tcfg.xmm*-1;
                ecfg.xmm=ecfg.xmm*-1;
            end
            if obj.planRelative(3)==2
                tcfg.ymm=tcfg.ymm*-1;
                ecfg.ymm=ecfg.ymm*-1;
            end
            if obj.planRelative(4)==1
                tcfg.zmm=tcfg.zmm*-1;
                ecfg.zmm=ecfg.zmm*-1;
            end
        end

        if fromspace==1 && tospace==3 % ACPC 2 MNI
            twarped=ea_acpc2mni(tcfg,{[obj.options.root,obj.options.patientname,filesep]});
            ewarped=ea_acpc2mni(ecfg,{[obj.options.root,obj.options.patientname,filesep]});
            t.target=twarped.WarpedPointMNI;
            t.entry=ewarped.WarpedPointMNI;
        elseif fromspace==1 && tospace==2 % ACPC 2 Native
            twarped=ea_acpc2mni(tcfg,{[obj.options.root,obj.options.patientname,filesep]});
            ewarped=ea_acpc2mni(ecfg,{[obj.options.root,obj.options.patientname,filesep]});
            t.target=twarped.WarpedPointNative;
            t.entry=ewarped.WarpedPointNative;
        elseif fromspace==3 && tospace==1 % MNI 2 ACPC
            twarped=ea_mni2acpc(tcfg,{[obj.options.root,obj.options.patientname,filesep]});
            ewarped=ea_mni2acpc(ecfg,{[obj.options.root,obj.options.patientname,filesep]});

            t.target=twarped.WarpedPointACPC;
            t.entry=ewarped.WarpedPointACPC;

        elseif fromspace==2 && tospace==1 % Native 2 ACPC
            twarped=ea_native2acpc(tcfg,{[obj.options.root,obj.options.patientname,filesep]});
            ewarped=ea_native2acpc(ecfg,{[obj.options.root,obj.options.patientname,filesep]});
            t.target=twarped.WarpedPointACPC;
            t.entry=ewarped.WarpedPointACPC;
        elseif fromspace==2 && tospace==3 % Native 2 MNI
            coords=[obj.target.target;obj.target.entry]';
            [obj.options,anats]=ea_assignpretra(obj.options);
            V=ea_open_vol([obj.options.root,obj.options.patientname,filesep,anats{1}]);
            coords=V.mat\[coords;1,1]; % go to voxel space in nativespace
            coords=ea_map_coords(coords, ...
                [obj.options.root,obj.options.patientname,filesep,anats{1}], ...
                [obj.options.root,obj.options.patientname,filesep,'inverseTransform'], ...
                '');
            t.target=coords(:,1)';
            t.entry=coords(:,2)';
        elseif fromspace==3 && tospace==2 % MNI 2 Native
            coords=[obj.target.target;obj.target.entry]';
            [obj.options,anats]=ea_assignpretra(obj.options);
            V=ea_open_vol([ea_space,obj.options.primarytemplate]);
            coords=V.mat\[coords;1,1]; % go to voxel space in MNI template
            src=[obj.options.root,obj.options.patientname,filesep,anats{1}]; % assign src image as primary anat image here.
            coords=ea_map_coords(coords, [ea_space,obj.options.primarytemplate], [obj.options.root,obj.options.patientname,filesep,'forwardTransform'],...
                src);

            t.target=coords(:,1)';
            t.entry=coords(:,2)';
        else % no change
            ea_busyaction('off',obj.controlH,'trajectory');
            return
        end

        if tospace==1 % to ACPC
            if obj.planRelative(2)==2
                t.target(1)=t.target(1)*-1;
                t.entry(1)=t.entry(1)*-1;
            end
            if obj.planRelative(3)==2
                t.target(2)=t.target(2)*-1;
                t.entry(2)=t.entry(2)*-1;
            end
            if obj.planRelative(4)==1
                t.target(3)=t.target(3)*-1;
                t.entry(3)=t.entry(3)*-1;
            end
        end

        obj.target=t;
        ea_synctrajectoryhandles(getappdata(obj.controlH,'chandles'),obj); % sync back control figure
        ea_save_electrode(obj);
        return
    end

    if ismember(evtnm,{'all','target','reco','hasPlanning','showMicro','relateMicro','planningAppearance','plan2elstruct_model','electrodeRelativeToPlan','color'}) % need to redraw planning fiducials:
        % planning fiducial
        if obj.showPlanning
            coords=ea_convertfiducials(obj,[obj.target.target;obj.target.entry]);
            tgt=coords(1,:); ent=coords(2,:);
            for dim=1:3
                traj(:,dim)=linspace(ent(dim),tgt(dim),10);
            end
            delete(obj.patchPlanning);

            % estimate pseudo-reconstruction (plan2elstruct):

            options=obj.options;
            options.elmodel=obj.plan2elstruct_model;
            options=ea_resolve_elspec(options);
            intraj=(ent-tgt)./norm(ent-tgt);

            if options.elspec.tipiscontact
                shift=-(obj.electrodeRelativeToPlan-1);
                markers.head=tgt+(((0+shift)*options.elspec.eldist)*intraj);
                markers.tail=tgt+(((3+shift)*options.elspec.eldist)*intraj);
            else
                if obj.electrodeRelativeToPlan>0
                    shift=-(obj.electrodeRelativeToPlan-1);
                    markers.head=tgt+(((0+shift)*options.elspec.eldist)*intraj);
                    markers.tail=tgt+(((3+shift)*options.elspec.eldist)*intraj);
                else % relative to tip
                    markers.head=tgt+(((options.elspec.tip_length/2)+(options.elspec.contact_length/2))*intraj);
                    markers.tail=tgt+(((3)*options.elspec.eldist)*intraj)+(((options.elspec.tip_length/2)+(options.elspec.contact_spacing))*intraj);
                end
            end

            [xunitv, yunitv] = ea_calcxy(markers.head, markers.tail);
            markers.x = markers.head + xunitv*(options.elspec.lead_diameter/2);
            markers.y = markers.head + yunitv*(options.elspec.lead_diameter/2);
            [coords_mm,trajectory,markers]=ea_resolvecoords(markers,options);
            obj.plan2elstruct(1).coords_mm=coords_mm;
            obj.plan2elstruct(1).trajectory=trajectory;
            obj.plan2elstruct(1).name='';
            obj.plan2elstruct(1).markers=markers;


            switch obj.planningAppearance
                case 'line'
                    [obj.patchPlanning, fv] = ea_plot3t(traj(:,1),traj(:,2),traj(:,3),obj.radius,obj.color,12,1);
                case 'electrode'
                    options=ea_defaultoptions(options);
                    options.sides=1;
                    options.colorMacroContacts=[];
                    obj.patchPlanning=ea_showelectrode(obj,'plan',options);
            end
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

%     if ismember(evtnm,{'color'}) % simply change color of patch
%         obj.patchPlanning.FaceVertexCData=repmat(obj.color,size(obj.patchPlanning.FaceVertexCData,1),1);
%     end

    if ismember(evtnm,{'showPlanning'}) && obj.hasPlanning
        for po=1:length({obj.patchPlanning.Visible})
        obj.patchPlanning(po).Visible=ea_bool2onoff(obj.showPlanning);
        end
    end

    if ismember(evtnm,{'all','elmodel','colorMacroContacts','elstruct'})
        if obj.showMacro
            try
                delete(obj.elpatch);
                delete(obj.ellabel);
            end
            poptions=obj.options;
            poptions.elmodel=obj.elmodel;
            obj.elstruct.elmodel=obj.elmodel;
            poptions.sides=obj.side;
            poptions.colorMacroContacts=obj.colorMacroContacts;
            el_render=getappdata(obj.plotFigureH,'el_render');
            if strcmp(evtnm,'elstruct')
                poptions.nowrite=1; % prevent from writing reconstruction to disk.
            end
            [obj.elpatch,obj.ellabel,obj.eltype]=ea_showelectrode(obj,'dbs',poptions);
            if isempty(el_render)
                clear el_render
            end
            el_render(obj.side).elpatch=obj.elpatch;
            setappdata(obj.plotFigureH,'el_render',el_render);
            if ~isempty(obj.ellabel)
                set(obj.ellabel,'Visible','off');
            end
        end
    end

    if ismember(evtnm,{'showMacro'}) && obj.hasMacro
        ea_elvisible([],[],obj.elpatch,1,obj.side,ea_bool2onoff(obj.showMacro),obj.options);
    end

    % add toggle button:
    if isfield(obj.elstruct, 'name') && ~isempty(obj.elstruct.name)
        ptname = obj.elstruct.name;
    else
        [~, ptname] = fileparts(fileparts(obj.options.root));
    end

    % Side label
    switch obj.side
        case 1
            elToggleSideLabel = 'Right';
        case 2
            elToggleSideLabel = 'Left';
        otherwise
            elToggleSideLabel = num2str(obj.side);
    end

    % Tooltip and icon
    if strcmp(obj.relateMicro, 'macro')
        elToggleTooltip = [ptname,' (',elToggleSideLabel,')'];
        elToogleIcon = ea_get_icn('electrode');
    else
        elToggleTooltip = [ptname,' (planning)'];
        elToogleIcon = ea_get_icn('electrode_planning');
    end

    % Add tag
    if strcmp(obj.options.leadprod, 'group') && isfield(obj.elstruct,'group')
        elToggleTag = ['Group: ', num2str(obj.elstruct.group), ...
               ', Patient: ', ptname, ', Side: ', elToggleSideLabel];
    elseif strcmp(obj.relateMicro, 'macro')
        elToggleTag = ['Patient: ', ptname, ', Side: ', elToggleSideLabel];
    elseif strcmp(obj.relateMicro, 'planning')
        elToggleTag = ['Patient: ', ptname, ', Planning'];
    end

    try ea_save_electrode(obj); end

    set(obj.toggleH, {'Parent','CData','TooltipString','Tag','OnCallback','OffCallback','State'},...
        {obj.htH, ...
        elToogleIcon, elToggleTooltip, elToggleTag,...
        {@ea_trajvisible,'on',obj}, {@ea_trajvisible,'off',obj}, ...
        ea_bool2onoff(any([obj.showPlanning,obj.showMacro,obj.showMicro]))});

    % Put all electrode toggles together
    isEleToggle = arrayfun(@(obj) ~isempty(regexp(obj.Tag, 'Patient', 'once')), allchild(obj.htH));
    EleToggleInd = find(isEleToggle)';
    if sum(diff(isEleToggle))~=0 % need to reorder
        notEleToggleNode = find(diff(isEleToggle))'+[1 0 1];
        newInd = [notEleToggleNode(1):notEleToggleNode(2),EleToggleInd,notEleToggleNode(3):length(isEleToggle)];
        obj.htH.Children = obj.htH.Children(newInd);
    end

    ea_busyaction('off',obj.controlH,'trajectory');
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
                            [obj.options.root,obj.options.patientname,filesep,'inverseTransform'], ...
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
                            [obj.options.root,obj.options.patientname,filesep,'forwardTransform'], ...
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


function fn=stripext(fn)
    [~,fn]=fileparts(fn);
end


function rightcallback(src, evt, obj)
    if evt.getButton() == 3
        ea_editfiducial(src,evt,obj)
    end
end


function ea_editfiducial(Hobj, evt, obj)
    obj.controlH=ea_trajectorycontrol(obj);
end


function ea_trajvisible(src, evt, onoff, obj)
    if getappdata(obj.plotFigureH,'altpressed') % hide all
        el_render=getappdata(obj.plotFigureH,'el_render');

        for el=1:length(el_render)
            el_render(el).showPlanning=ea_bool2onoff('off');
            el_render(el).showMacro=ea_bool2onoff(onoff);
            el_render(el).showMicro=ea_bool2onoff('off');
        end
    elseif getappdata(obj.plotFigureH,'cmdpressed') % rightclick
        ea_editfiducial(src,evt,obj)
    else
        if ~isempty(obj.controlH) ...
           && ~strcmp(evt.Source.Type, 'uitoggletool') % From control figure
            chandles = getappdata(obj.controlH,'chandles');
            obj.showPlanning = get(chandles.showPlanning,'Value');
            obj.showMacro = get(chandles.showMacro,'Value');
        else % From toogle button
            switch onoff
                case 'off'
                    if strfind(evt.Source.Tag, 'Planning')
                        obj.showPlanning=0;
                    else
                        obj.showMacro=0;
                        obj.showMicro=0;
                    end

                    if ~isempty(obj.controlH)
                        chandles = getappdata(obj.controlH,'chandles');
                        if strfind(evt.Source.Tag, 'Planning')
                            set(chandles.showPlanning, 'Value', 0);
                        else
                            set(chandles.showMacro, 'Value', 0);
                            set(chandles.showMicro, 'Value', 0);
                        end
                    end
                case 'on'
                    if strfind(evt.Source.Tag, 'Planning')
                        obj.showPlanning=1;
                    else
                        obj.showMacro=1;
                    end

                    if ~isempty(obj.controlH)
                        chandles = getappdata(obj.controlH,'chandles');
                        if strfind(evt.Source.Tag, 'Planning')
                            set(chandles.showPlanning, 'Value', 1);
                        else
                            set(chandles.showMacro, 'Value', 1);
                        end
                    end
            end
        end
    end
end
