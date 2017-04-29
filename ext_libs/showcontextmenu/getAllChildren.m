function children = getAllChildren( hObject )
%GETALLCHILDREN 
%this returns all children as if get(0,'ShowHiddenHandles') == 'on'

old=get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');

children = get(hObject,'children');

set(0,'ShowHiddenHandles',old);

end

