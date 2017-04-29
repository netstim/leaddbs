function pos = getCurrentPoint( hFigure,units )
%GETCURRENTPOINT returns get(hFigure,'currentpoint') with proper set(hObject,'units',units);

old_units=get(hFigure,'units');
set(hFigure,'units',units);
pos=get(hFigure,'currentpoint');
set(hFigure,'units',old_units);

end

