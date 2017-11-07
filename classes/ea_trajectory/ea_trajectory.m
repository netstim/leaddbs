classdef ea_trajectory < handle
    % Trajectory class to handle trajectories in lead dbs resultfig / 3D Matlab figures
    % A. Horn
    
    properties (SetObservable)
        reco % reconstruction struct as saved in ea_reconstruction.mat files, available in native, MNI and AC/PC spaces
        target % target and entrypoints as used in surgical planning
        alpha=0.7 % alpha of Planning patch
        radius=0.2 % radius of Planning line
        color=[0.8,0.3,0.2] % color of Planning patch
        options % lead-dbs options struct
        planRelative=[2,1,1,1,1] % First entry: AC=1, MCP=2, PC=3; Second entry: Right=1, Left=2; Third entry: Anterior=1, Posterior=2; Fourth entry: Ventral=1; Dorsal=2; Last entry: ACPC=1, native=2, MNI/Template=3
        hasPlanning % determines if object has information to show a fiducial
        hasMacro % determines if object has information to show a macroelectrode
        relateMicro='macro' % determines if microelectrodes shown should be related to planning Fiducial ('planning') or Macroelectrodes ('macro')
        showPlanning=1 % show planning fiducial
        showMacro=0 % show definitive DBS / macro electrode
        showMicro=0 % show microelectrodes
        controlH % handle to trajectory control figure
        plotFigureH % handle of figure on which to plot
        patchMacro % handle of macroelectrode patch
        patchPlanning % handle of planning fiducial patch
        patchMicro % handle of microelectrodes
        toggleH % togglebutton handle that will open planning fiducial control
        htH % handle for toggle toolbar on which toggleH is displayed
    end
    
    methods
        function obj=ea_trajectory(pobj) % generator function
            try
                obj.plotFigureH=pobj.plotFigureH;
            catch
                obj.plotFigureH=gcf;
            end
            
            
            obj.htH=getappdata(obj.plotFigureH,'addht');
            
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
                obj.reco=pobj.reco;
            catch
                obj.reco=struct;
            end
            if ~exist('pobj','var') % create blank trajectory with planning fiducial only
                obj.hasPlanning=1;
                obj.hasMacro=0;
            else % determine if fiducial and macro information is available
                obj.hasMacro=~isempty(obj.reco);
                obj.hasPlanning=~isempty(obj.target);
            end
            %% initialize further content fields based on given struct if given or else empty / random vars
            if obj.hasPlanning
                try % target
                    obj.target=pobj.target;
                catch
                    obj.target.entry=[20,20,50];
                    obj.target.target=[12.02,-1.53,1.91];
                    obj.target.offset=0;
                    obj.target.hemisphere=1; % right
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
            
            
            obj.options=getappdata(obj.plotFigureH,'options');
            if isempty(obj.options)
                ea_warning('Patient information not available.');
            end
            obj.toggleH=uitoggletool;
            
            % Get the underlying java object using findobj
            jtoggle = findjobj(obj.toggleH);
            
            % Specify a callback to be triggered on any mouse release event
            set(jtoggle, 'MouseReleasedCallback', {@rightcallback,obj})
            update_trajectory(obj);
            addlistener(obj,'showPlanning','PostSet',...
                @ea_trajectory.changeevent);
            addlistener(obj,'showMacro','PostSet',...
                @ea_trajectory.changeevent);
            addlistener(obj,'showMicro','PostSet',...
                @ea_trajectory.changeevent);
            addlistener(obj,'color','PostSet',...
                @ea_trajectory.changeevent);
            
            addlistener(obj,'target','PostSet',...
                @ea_trajectory.changeevent);
            
            addlistener(obj,'alpha','PostSet',...
                @ea_trajectory.changeevent);
            
            addlistener(obj,'target','PostSet',...
                @ea_trajectory.changeevent);
            
            if (exist('pobj','var') && isfield(pobj,'openedit') && pobj.openedit) || ~exist('pobj','var')
                ea_trajectorycontrol(obj)
            end
            
        end
        
        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
        
        function obj=update_trajectory(obj,evtnm) % update ROI
            if ~exist('evtnm','var')
                evtnm='all';
            end
            set(0,'CurrentFigure',obj.plotFigureH);
            if ismember(evtnm,{'all','target','reco','color'}) && obj.hasPlanning % need to redraw planning fiducials:
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
            end
            if ismember(evtnm,{'showPlanning'}) && obj.hasPlanning
                                obj.patchPlanning.Visible=ea_bool2onoff(obj.showPlanning);
            end
            % add toggle button:
            set(obj.toggleH,...
                {'Parent','CData','TooltipString','OnCallback','OffCallback','State'},...
                {obj.htH,ea_get_icn('atlas',obj.color),'Trajectory',{@ea_trajvisible,'on',obj},{@ea_trajvisible,'off',obj},ea_bool2onoff(obj.showPlanning)});
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
                        if obj.planRelative(3)==2; cfg.ymm=-cfg.xmm; end
                        if obj.planRelative(4)==2; cfg.zmm=-cfg.xmm; end
                        switch obj.options.native
                            case 1 % need to convert from AC/PC to native
                                
                            case 0 % need to convert from AC/PC to template
                                ccoords(coord,:)=ea_acpc2mni(cfg,{[obj.options.root,obj.options.patientname,filesep]});
                        end
                    case 2 % planning in native space
                        switch obj.options.native
                            case 1 % leave coords as they are
                                
                            case 0 % need to convert from native to MNI
                                
                        end
                    case 3 % planning in template space
                        switch obj.options.native
                            case 1 % need to convert from MNI to native
                                
                            case 0 % leave coords as they are
                                
                        end
                        
                end
            end
        end
        
        function rightcallback(src, evnt,obj)
            if evnt.getButton() == 3
                ea_editfiducial(src,evnt,obj)
            end
        end
        
        function ea_editfiducial(Hobj,evt,obj)
            obj.controlH=ea_trajectorycontrol(obj);
            
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
    end
    
end

function ea_trajvisible(~,~,onoff,obj)
switch onoff
    case 'on'
        obj.showMacro=1;
        obj.showMicro=1;
        obj.showPlanning=1;
    case 'off'
        obj.showMacro=0;
        obj.showMicro=0;
        obj.showPlanning=0;
end
end
