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

% Update anatomy control GUI
awin = getappdata(resultfig, 'awin');
handle = guidata(awin);
handle.xval.String = num2str(xx);
handle.yval.String = num2str(yy);
handle.zval.String = num2str(zz);

% Update slices
ea_anatomyslices(resultfig, togglestates, options, []);
