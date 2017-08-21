function ea_setaxisorigin(options, resultfig)
% attempt to zoom at the center of the electrode reconstruction.
% This is a highly optional step and only works if a patient is
% selected, thus whole step is in try/end brackets for now.

% calculate the center of the contacts
coords = ea_load_reconstruction(options);
origin = mean(cell2mat(cellfun(@mean, coords, 'Uniformoutput',0)'), 1);

% set the origin (zoom point) at the center of the contacts
xrange = max(abs(resultfig.CurrentAxes.XLim - origin(1)));
yrange = max(abs(resultfig.CurrentAxes.YLim - origin(2)));
zrange = max(abs(resultfig.CurrentAxes.ZLim - origin(3)));

axis([origin(1) + [-1,1] * xrange, origin(2) + [-1,1] * yrange, origin(3) + [-1,1] * zrange]);

% zoon out
zoom(3)
