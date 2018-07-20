function ea_setbackdrop(backdrop,options)

H = findall(0,'type','figure');
resultfig = H(~cellfun(@isempty,strfind({H(:).Name},{'Electrode-Scene'})));
togglestates=getappdata(resultfig,'togglestates');
togglestates.template=backdrop;
if ~exist('options','var')
    options=struct;
end
setappdata(resultfig,'togglestates',togglestates);

ea_anatomyslices(gcf,...
    togglestates,...
    options,[]);