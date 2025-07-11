function dndobj = ea_bind_dragndrop(target, DropFileFcn, DropStringFcn)
% Bind MATLAB uicontrol or figure with drag-drop event

if isMatlabVer('>=', [25,1])
    % Won't work since JavaFrame has been removed since R2025a
    return;
end

% check if target is a whole figure
if isa(target,'matlab.ui.Figure')
    target = ea_getJavaContentPane(target);
end

% initialize the dndcontrol class
dndcontrol.initJava();
dndobj = dndcontrol(target);

% optionally set callback function
if nargin >=2
    dndobj.DropFileFcn = DropFileFcn;
end

if nargin >=3
    dndobj.DropStringFcn = DropStringFcn;
end
