function ea_mercontrol_updatemarkers(handles)
% Replace the string value in handles.popupmermarkers_right _left
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
mermarkers = getappdata(resultfig, 'mermarkers');

popup_str_array = cell(1, length(mermarkers));
for marker_ix = 1:length(mermarkers)
    marker = mermarkers(marker_ix);
    popup_str_array{marker_ix} = sprintf('%d. %s', marker_ix, marker.tag.string);
end
side_str_array = {mermarkers.side};

for side_str = {'right', 'left'}
    this_str_array = cat(2, {'none selected...'},...
        popup_str_array(strcmpi(side_str_array, side_str{:})));
    set(handles.(['popupmermarkers_' side_str{:}]),...
        'Visible', 'off', 'Value', 1,...
        'String', this_str_array);
end

if isempty(mermarkers)
    enab_str = 'off';
else
    enab_str = 'on';
end
handles.undomarker.Enable = enab_str;
handles.togglemarkertags.Enable = enab_str;
handles.exportmarkers.Enable = enab_str;

markerhistory = getappdata(resultfig, 'markerhistory');
if isempty(markerhistory)
    handles.redomarker.Enable = 'off';
else
    handles.redomarker.Enable = 'on';
end

if strcmpi(merstruct.tag.visible, 'off')
    set(handles.togglemarkertags, 'Value', 0);
else
    set(handles.togglemarkertags, 'Value', 1);
end