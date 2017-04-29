function position = getPositionOnFigure( hObject,units )
%GETPOSITIONONFIGURE returns absolute position of object on a figure
%
% (since get(hObject,'Position') returns position relative to hObject's
%  parent)

hObject_pos=getRelPosition(hObject,units);


parent = get(hObject,'Parent');
parent_type = get(parent,'Type');

if isequal(parent_type,'figure')
  position = hObject_pos;
  return;

else
    parent_pos = getPositionOnFigure( parent,units );
    position = relativePos2absolutePos(hObject_pos,parent_pos,units);    
end

