function v=ea_view(varargin);
% query and store views
if nargin
    v=varargin{1};
    view([v.az,v.el]);
    camva(v.camva);
    camup(v.camup);
    camproj(v.camproj);
    camtarget(v.camtarget);
    campos(v.campos); 
else
    [v.az,v.el]=view;
    v.camva=camva;
    v.camup=camup;
    v.camproj=camproj;
    v.camtarget=camtarget;
    v.campos=campos;  
end