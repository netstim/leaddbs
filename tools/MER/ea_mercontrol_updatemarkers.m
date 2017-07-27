function ea_mercontrol_updatemarkers(handles)
% Replace the string value in handles.popupmermarkers_right _left
merstruct = getappdata(handles.mercontrolfig, 'merstruct');

popup_str_array = cell(1, length(merstruct.Markers));
for marker_ix = 1:length(merstruct.Markers)
    marker = merstruct.Markers(marker_ix);
    tag_string = sprintf('%s.%s.%s: %.2fmm',...
            marker.side(1), marker.tract_label(1:3), marker.type(1:3), marker.depth);
    popup_str_array{marker_ix} = sprintf('%d. %s', marker_ix, tag_string);
end
side_str_array = {merstruct.Markers.side};

for side_str = {'right', 'left'}
    this_str_array = cat(2, {'none selected...'},...
        popup_str_array(strcmpi(side_str_array, side_str{:})));
    set(handles.(['popupmermarkers_' side_str{:}]),...
        'Visible', 'on', 'Value', 1,...
        'String', this_str_array);
end

if isempty(merstruct.Markers)
    enab_str = 'off';
else
    enab_str = 'on';
end
handles.undomarker.Enable = enab_str;
handles.togglemarkertags.Enable = enab_str;

if isempty(merstruct.MarkersHistory)
    handles.redomarker.Enable = 'off';
else
    handles.redomarker.Enable = 'on';
end