function ea_setplanes(xx,yy,zz,options)

H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
togglestates=getappdata(resultfig,'togglestates');
setXYZ=[xx,yy,zz];
togglestates.xyztoggles=~isnan(setXYZ);
setXYZ(~togglestates.xyztoggles)=togglestates.xyzmm(~togglestates.xyztoggles);
togglestates.xyzmm=setXYZ;
togglestates.refreshview=1;
if ~exist('options','var')
    options=struct;
end
ea_anatomyslices(gcf,...
    togglestates,...
    options,[]);