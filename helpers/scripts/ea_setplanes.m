function ea_setplanes(xx,yy,zz,options,template)

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
if exist('template','var')
    togglestates.template=template;
    togglestates.refreshcuts=1;
end

% Set togglestates
setappdata(resultfig,'togglestates',togglestates);

% Update anatomy control GUI
awin = getappdata(resultfig, 'awin');
if isvalid(awin) % check if not deleted (closed)
    handle = guidata(awin);
    handle.xval.String = num2str(xx);
    handle.yval.String = num2str(yy);
    handle.zval.String = num2str(zz);
end

togglestates.xyztoggles=~isnan(setXYZ);

% Update slices
ea_anatomyslices(resultfig, togglestates, options, []);
