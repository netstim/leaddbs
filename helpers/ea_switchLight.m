function ea_switchLight(figHandle,viz)
% Switch on/off the light in 3D visualization window

if nargin<2
    viz=1;
end

if numel(viz) == 1
    viz = repmat(viz, 1, 4);
end

set(0,'CurrentFigure',figHandle);

CamLight = getappdata(figHandle, 'CamLight');
CeilingLight = getappdata(figHandle, 'CeilingLight');
RightLight = getappdata(figHandle, 'RightLight');
LeftLight = getappdata(figHandle, 'LeftLight');

if ~isempty(CamLight)
    if viz(1) == 1
        set(CamLight,'Visible','on');
    else
        set(CamLight,'Visible','off');
    end
end

if ~isempty(CeilingLight)
    if viz(2) == 1
        set(CeilingLight,'Visible','on');
    else
        set(CeilingLight,'Visible','off');
    end
end

if ~isempty(RightLight)
    if viz(3) == 1
        set(RightLight,'Visible','on');
    else
        set(RightLight,'Visible','off');
    end
end

if ~isempty(LeftLight)
    if viz(4) == 1
        set(LeftLight,'Visible','on');
    else
        set(LeftLight,'Visible','off');
    end
end
