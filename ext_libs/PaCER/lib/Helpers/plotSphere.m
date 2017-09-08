%% plotSphere - plot a sphere for given point
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function graphicsHandle = plotSphere(point, diameter, color, parent, varargin)
    if(nargin < 5)
        varargin = {};
    end
    if(nargin < 4)
        parent = gca;
    end
    if(nargin < 3)
        color = 'y';
    end

    [x,y,z] = sphere(); % 10x10 faces (default: 20x20 facess)
    
    x = x .* (diameter/2)+ point(1);
    y = y .* (diameter/2)+ point(2);
    z = z .* (diameter/2)+ point(3);
    
    graphicsHandle = surf(x,y,z, 'FaceColor', color, 'EdgeColor', 'none', 'Parent', parent, varargin{:});
    daspect([1 1 1]);
    lighting gouraud;
    material shiny;
end