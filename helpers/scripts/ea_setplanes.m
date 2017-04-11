function ea_setplanes(xx,yy,zz,options)

H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
togglestates=getappdata(resultfig,'togglestates');
togglestates.xyzmm=[xx,yy,zz];
togglestates.refreshview=1;
if ~exist('options','var')
    options=struct;
end
ea_anatomyslices(gcf,...
    togglestates,...
    options,[]);