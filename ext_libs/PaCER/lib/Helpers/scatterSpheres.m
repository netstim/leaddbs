%% scatterSpheres - plot scattered points as spheres
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function graphicsHandle = scatterSpheres(points, diameter, color, varargin)
    if(nargin < 4)
        varargin = {};
    end
    if(nargin < 3)
        color = 'y';
    end
    
    nm = size(points);
    if(nm(1) < nm(2))
        points = points';
    end
    graphicsHandle = hggroup;
    for i = 1:length(points);
        plotSphere(points(i,:), diameter, color, graphicsHandle, varargin{:});
    end
end