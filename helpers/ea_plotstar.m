function graphicsHandle = ea_plotstar(point, diameter, color, edgecolor, alpha)

    if ~exist('edgecolor','var')
        edgecolor = 'none';
    end
    if ~exist('alpha','var')
        alpha = 1;
    end

    [x,y,z] = sphere(30);

    % Convert to spherical coordinates
    [az,el,r] = cart2sph(x,y,z);
    
    nAz = 4;      % number of spikes around
    nEl = 4;      % number of spikes vertically
    amp = 0.8;   % spikiness strength

    spikes = 1 + amp * cos(nAz*az) .* cos(nEl*el);
    % Create spikes
    %spikes = 1 + 0.35 * cos(5*az) .* cos(5*el);
    %spikes = 1 + 0.20 * cos(5*az) .* cos(5*el);
    %spikes = 1 + 0.22 * cos(3*az) .* cos(3*el);
    r = r .* spikes;

    % Back to cartesian
    [x,y,z] = sph2cart(az, el, r);

    % Scale + translate
    x = x * (diameter/2) + point(1);
    y = y * (diameter/2) + point(2);
    z = z * (diameter/2) + point(3);

    graphicsHandle = surf(x,y,z, ...
        'FaceColor', color, ...
        'EdgeColor', edgecolor, ...
        'FaceAlpha', alpha);

    daspect([1 1 1])
    lighting gouraud
    material shiny
end