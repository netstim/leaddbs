function ea_setplanes(xx,yy,zz,options)

H = findall(0,'type','figure');
resultfig = H(contains({H(:).Name},{'Electrode-Scene'}));
resultfig = resultfig(1); % take the first if there are many.
togglestates=getappdata(resultfig,'togglestates');
setXYZ=[xx,yy,zz];
togglestates.xyztoggles=~isnan(setXYZ);
setXYZ(~togglestates.xyztoggles)=togglestates.xyzmm(~togglestates.xyztoggles);
togglestates.xyzmm=setXYZ;
togglestates.refreshview=1;
if ~exist('options','var')
    options=struct;
end
setappdata(resultfig,'togglestates',togglestates);

ea_anatomyslices(resultfig, togglestates, options, []);
