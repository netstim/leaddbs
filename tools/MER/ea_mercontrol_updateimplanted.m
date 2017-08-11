function ea_mercontrol_updateimplanted(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, side_ids, ~] = ea_detsidestr(side_str);
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
ui_tags = {handles.mainuipanel.Children.Tag};
for sid = side_ids
    side_str = side_strs{sid};
    b_item = strcmpi(handles.(['popupimplantedtract_' side_str]).String,...
        merstruct.DBSImplants(sid).implanted_tract_label);
    set(handles.(['popupimplantedtract_' side_str]),...
        'Value', find(b_item));
    set(handles.(['editimplanteddepth_' side_str]),...
        'String', num2str(merstruct.DBSImplants(sid).depth));
%         set(handles.(['popupmermarkers_' side_str]),...
%             'Visible', 'off', 'String', '', 'Value', 1);
    % TODO: Recording elements handles.(['popupmedial_' side_str])
end
for traj_ix = 1:length(merstruct.MERTrajectories)
    traj = merstruct.MERTrajectories(traj_ix);
    set(handles.mainuipanel.Children(strcmpi(ui_tags,...
        ['pos', traj.label, '_', traj.side])), 'String', num2str(traj.depth));
end