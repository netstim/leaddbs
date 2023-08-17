%% plotSphere - plot a sphere for given point
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu

% adapted for use in lead-dbs

function graphicsHandle = ea_plotsphere(point, diameter, color, edgecolor, alpha)

    if ~exist('edgecolor','var')
        edgecolor='none';
    end

    if ~exist('alpha','var')
        alpha=1;
    end

    [x,y,z] = sphere(10); % 10x10 faces (default: 20x20 facess)
    
    x = x .* (diameter/2)+ point(1);
    y = y .* (diameter/2)+ point(2);
    z = z .* (diameter/2)+ point(3);
    
    graphicsHandle = surf(x,y,z, 'FaceColor', color, 'EdgeColor', edgecolor, 'FaceAlpha', alpha);
    daspect([1 1 1]);
    lighting gouraud;
    material shiny;
end