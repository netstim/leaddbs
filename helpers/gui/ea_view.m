function v=ea_view(varargin)
% query and store views
if nargin
    v=varargin{1};
    view([v.az,v.el]);
    camva(v.camva);
    camup(v.camup);
    camproj(v.camproj);
    camtarget(v.camtarget);
    campos(v.campos);
    if nargin>1
        resultfig=varargin{2};
    else
        resultfig=gcf;
    end
else
    resultfig=gcf;
    [v.az,v.el]=view;
    v.camva=camva;
    v.camup=camup;
    v.camproj=camproj;
    v.camtarget=camtarget;
    v.campos=campos;
end

set(0,'CurrentFigure',resultfig)
allAxes = findall(resultfig,'type','axes');
set(resultfig,'CurrentAxes',allAxes(1));
ea_zoomcenter(resultfig.CurrentAxes, v.camtarget);
