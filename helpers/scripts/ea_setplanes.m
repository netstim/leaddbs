function ea_setplanes(xx,yy,zz,options)

% Get 3D visualization figure handle
H = findall(0,'type','figure');
resultfig = H(contains({H(:).Name},{'Electrode-Scene'}));
resultfig = resultfig(1); % take the first if there are many.

% Set non-nan values to xyzmm
togglestates = getappdata(resultfig,'togglestates');
setXYZ = [xx,yy,zz];
togglestates.xyzmm(~isnan(setXYZ)) = setXYZ(~isnan(setXYZ));

% Flag to refresh view
togglestates.refreshview=1;

% Options for ea_anatomyslices input
if ~exist('options','var')
    options=struct;
end

% Set togglestates
setappdata(resultfig,'togglestates',togglestates);


% Update slices
ea_anatomyslices(resultfig, togglestates, options, []);
