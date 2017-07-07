function ea_mercontrol_updateimplanted(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, side_ids, ~] = ea_detsidestr(side_str);
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
ui_tags = {handles.mainuipanel.Children.Tag};
for sid = side_ids
    side_str = side_strs{sid};
    set(handles.(['popupimplantedtract_' side_str]),...
        'Value', merstruct.implant_idx(sid) + 1);
    set(handles.(['editimplanteddepth_' side_str]),...
        'String', num2str(merstruct.implant_depth(sid)));
%         set(handles.(['popupmermarkers_' side_str]),...
%             'Visible', 'off', 'String', '', 'Value', 1);
    % TODO: Recording elements handles.(['popupmedial_' side_str])
    for tract_info = merstruct.tract_info
        if isfield(merstruct, 'currentmer')
            val = merstruct.currentmer.(tract_info.label).dist(sid);
        else
            val = 0;
        end
        set(handles.mainuipanel.Children(strcmpi(ui_tags,...
            ['pos' tract_info.label '_' side_str])),...
            'String', num2str(val));
    end
end