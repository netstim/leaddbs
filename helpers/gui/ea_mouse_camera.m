function [] = ea_mouse_camera(hfig)
% Execute this function to a Figure hfig in order to set camera movements
% to the mouse.
%
% Left click: rotate
% Wheel click / Left click + shift: move
% Right click / Left click + ctrl: zoom


figLastPoint = []; % global variable to store previous cursor position
zoomFactor = 10;
panFactor = 2;
orbitFactor = 7;

set(hfig, 'WindowButtonDownFcn', @down_fcn);
set(hfig, 'WindowButtonUpFcn', @up_fcn);
set(hfig, 'WindowScrollWheelFcn', {@(src, evt) zoom_fcn(evt)});

    function [] = zoom_fcn(evt)
        if evt.VerticalScrollCount > 0
            camzoom(1 + evt.VerticalScrollCount / zoomFactor)
        else
            camzoom(1 / (1 + abs(evt.VerticalScrollCount) / zoomFactor))
        end
    end


    function [] = down_fcn(hfig, evt)
        clickType = evt.Source.SelectionType;
        set(hfig, 'WindowButtonMotionFcn',{@(src, evt) motion_callback(src, clickType)});
        
        % set cursor type
        switch clickType
            case 'normal'
                setptr(gcf, 'rotate');
            case 'alt'
                setptr(gcf, 'hand');
            case 'extend'
                setptr(gcf, 'glass');
            case 'open'
                try
                    set_defaultview;
                end  
        end
    end

    function set_defaultview
        % get stored default view preferences and call ea_defaultview
        prefs = ea_prefs;
        v = prefs.machine.view;
        togglestates = prefs.machine.togglestates;
        ea_defaultview_transition(v,togglestates);
        ea_defaultview(v,togglestates);
    end

    function [] = motion_callback(hfig, clickType)
        % from matlab's CameraToolBarManager.m
        currpt = get(hfig,'CurrentPoint');
        try
            pt = matlab.graphics.interaction.internal.getPointInPixels(hfig,currpt(1:2));
        catch % old matlab version
            pt = currpt;
        end
        if isempty(figLastPoint)
            figLastPoint = pt;
        end
        deltaPix  = pt-figLastPoint;
        figLastPoint = pt;
        
        switch clickType
            case 'normal'
                orbitPangca(deltaPix/orbitFactor, 'o');
            case 'alt'
                dollygca(deltaPix/panFactor);
            case 'extend'
                zoomgca(deltaPix);
        end
        
    end

    function dollygca(xy)
        % from matlab's CameraToolBarManager.m
        haxes = gca;
        camdolly(haxes,-xy(1), -xy(2), 0, 'movetarget', 'pixels')
        drawnow
    end

    function orbitPangca(xy, mode)
        % from matlab's CameraToolBarManager.m
        %mode = 'o';  orbit
        %mode = 'p';  pan
        
        %coordsystem = lower(hObj.coordsys);
        coordsystem = 'z';
        
        haxes = gca;
        
        if coordsystem(1)=='n'
            coordsysval = 0;
        else
            coordsysval = coordsystem(1) - 'x' + 1;
        end
        
        xy = -xy;
        
        if mode=='p' % pan
            panxy = xy*camva(haxes)/500;
        end
        
        if coordsysval>0
            d = [0 0 0];
            d(coordsysval) = 1;
            
            up = camup(haxes);
            upsidedown = (up(coordsysval) < 0);
            if upsidedown
                xy(1) = -xy(1);
                d = -d;
            end
            
            % Check if the camera up vector is parallel with the view direction;
            % if not, set the up vector
            try
                check = any(matlab.graphics.internal.CameraToolBarManager.crossSimple(d,campos(haxes)-camtarget(haxes)));
            catch % Matlab 2017
                check = any(crossSimple(d,campos(haxes)-camtarget(haxes)));
            end
            if check
                camup(haxes,d)
            end
        end
        
        flag = 1;
        
        %while sum(abs(xy))> 0 && (flag || hObj.moving) && ishghandle(haxes)
        while sum(abs(xy))> 0 && (flag) && ishghandle(haxes)
            flag = 0;
            if ishghandle(haxes)
                if mode=='o' %orbit
                    if coordsysval==0 %unconstrained
                        camorbit(haxes,xy(1), xy(2), coordsystem)
                    else
                        camorbit(haxes,xy(1), xy(2), 'data', coordsystem)
                    end
                else %pan
                    if coordsysval==0 %unconstrained
                        campan(haxes,panxy(1), panxy(2), coordsystem)
                    else
                        campan(haxes,panxy(1), panxy(2), 'data', coordsystem)
                    end
                end
                %updateScenelightPosition(hObj,haxes);
                %localDrawnow(hObj);
                drawnow
            end
        end
    end

    function zoomgca(xy)
        % from matlab's CameraToolBarManager.m
        
        haxes = gca;
        
        q = max(-.9, min(.9, sum(xy)/70));
        q = 1+q;
        
        % heuristic avoids small view angles which will crash on Solaris
        MIN_VIEW_ANGLE = .001;
        MAX_VIEW_ANGLE = 75;
        vaOld = camva(gca);
        camzoom(haxes,q);
        va = camva(haxes);
        %If the act of zooming puts us at an extreme, back the zoom out
        if ~((q>1 || va<MAX_VIEW_ANGLE) && (va>MIN_VIEW_ANGLE))
            set(haxes,'CameraViewAngle',vaOld);
        end
        
        drawnow
    end

    function [] = up_fcn(hfig, evt)
        % Adjust light after rotating
        if strcmp(evt.Source.SelectionType, 'normal')
            ea_readjustlight(hfig);
        end
        % reset motion and cursor
        set(hfig,'WindowButtonMotionFcn',[]);
        figLastPoint = [];
        setptr(gcf, 'arrow');
    end

    function c=crossSimple(a,b)
        c(1) = b(3)*a(2) - b(2)*a(3);
        c(2) = b(1)*a(3) - b(3)*a(1);
        c(3) = b(2)*a(1) - b(1)*a(2);
    end

    function ea_readjustlight(hfig)
        % get light handles
        RightLight = getappdata(hfig, 'RightLight');
        LeftLight = getappdata(hfig, 'LeftLight');
        CeilingLight = getappdata(hfig, 'CeilingLight');
        CamLight = getappdata(hfig, 'CamLight');

        prefs=ea_prefs;

        try
            camlight(CamLight, 'headlight'); % move light object.
        end
        try
            set(CeilingLight, 'Position', [0 0 10], 'style', 'local', 'Color', prefs.d3.ceilinglightcolor); % not modifiable, infinite light.
        end
        try
            set(RightLight, 'Position', [-100 0 0], 'style', 'infinite', 'Color', prefs.d3.rightlightcolor); % not modifiable, infinite light.
        end
        try
            set(LeftLight, 'Position', [100 0 0], 'style', 'infinite', 'Color', prefs.d3.leftlightcolor); % not modifiable, infinite light.
        end
    end
end
