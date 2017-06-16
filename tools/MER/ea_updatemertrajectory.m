function ea_updatemertrajectory(handles,trajectory,dist,tag)
% ea_updatemertrajectory(handles,trajectory,dist,tag)
resultfig=getappdata(handles.mercontrolfig,'resultfig');
% Update position in resultfig
% XData = get(getappdata(resultfig,tag),'XData');
% YData = get(getappdata(resultfig,tag),'YData');
% ZData = get(getappdata(resultfig,tag),'ZData');
[~,side,track]=ea_detsidestr(tag);
h = getfield(getappdata(resultfig,'merhandles'),track);
h = h{side};
set(h,'XData',trajectory(:,1)')
set(h,'YData',trajectory(:,2)')
set(h,'ZData',trajectory(:,3)')
setappdata(resultfig,tag,h)
set(handles.(tag),'Value',1)
setappdata(resultfig,'keymer',tag)
newdiststr = num2str(str2double(get(handles.(['pos',tag(4:end)]),'String'))+dist);
set(handles.(['pos',tag(4:end)]),'String',newdiststr)