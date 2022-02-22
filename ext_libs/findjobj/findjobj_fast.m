function jControl = findjobj_fast(hControl, jContainer) %#ok<*JAVFM,*JAPIMATHWORKS>
    if double(get(hControl,'Parent'))~=0  % avoid below for figure handles
        try jControl = hControl.getTable; return, catch, end  % fast bail-out for old uitables
        try jControl = hControl.JavaFrame.getGUIDEView; return, catch, end  % bail-out for HG2 matlab.ui.container.Panel
    end
    oldWarn = warning;
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved');
    try oc = onCleanup(@()warning(oldWarn)); catch, end % Gracefully restore warnings on return (suggested by T. Carpenter 2021-08-17)
    if nargin < 2 || isempty(jContainer)
        % Use a HG2 matlab.ui.container.Panel jContainer if the control's parent is a uipanel
        try
            hParent = get(hControl,'Parent');
        catch
            % Probably indicates an invalid/deleted/empty handle
            jControl = [];
            return
        end
        try jContainer = hParent.JavaFrame.getGUIDEView; catch, jContainer = []; end
    end
    if isempty(jContainer)
        if isequal(double(hControl),0)  % root handle => return jDesktop
            try
                jControl = com.mathworks.mde.desk.MLDesktop.getInstance;
                if isempty(jControl), error('retry'), end
            catch
                jControl = com.mathworks.mlservices.MatlabDesktopServices.getDesktop;
            end
            return
        end
        hFig = ancestor(hControl,'figure');
        jf = get(hFig, 'JavaFrame');
        if isequal(hFig,hControl) % speedup suggested by T. Carpenter 2021-08-17
            jControl = jf;
            return
        end
        jContainer = jf.getFigurePanelContainer.getComponent(0);
    end
    warning(oldWarn);
    jControl = [];
    counter = 20;  % 2018-09-21 speedup (100 x 0.001 => 20 x 0.005) - Martin Lehmann suggestion on FEX 2016-06-07
    specialTooltipStr = '!@#$%^&*';
    try  % Fix for R2018b suggested by Eddie (FEX comment 2018-09-19)
        tooltipPropName = 'TooltipString';
        oldTooltip = get(hControl,tooltipPropName);
        set(hControl,tooltipPropName,specialTooltipStr);
    catch
        tooltipPropName = 'Tooltip';
        oldTooltip = get(hControl,tooltipPropName);
        set(hControl,tooltipPropName,specialTooltipStr);
    end
    while isempty(jControl) && counter>0
        counter = counter - 1;
        pause(0.005);
        jControl = findTooltipIn(jContainer, specialTooltipStr);
    end
    set(hControl,tooltipPropName,oldTooltip);
    try jControl.setToolTipText(oldTooltip); catch, end
    try jControl = jControl.getParent.getView.getParent.getParent; catch, end  % return JScrollPane if exists
end

function jControl = findTooltipIn(jContainer, specialTooltipStr)
    try
        jControl = [];  % Fix suggested by H. Koch 11/4/2017
        tooltipStr = jContainer.getToolTipText;
        %if strcmp(char(tooltipStr),specialTooltipStr)
        if ~isempty(tooltipStr) && tooltipStr.startsWith(specialTooltipStr)  % a bit faster
            jControl = jContainer;
        else
            for idx = 1 : jContainer.getComponentCount
                jControl = findTooltipIn(jContainer.getComponent(idx-1), specialTooltipStr);
                if ~isempty(jControl), return; end
            end
        end
    catch
        % ignore
    end
end
