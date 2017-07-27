function ea_resultfig_updatetrajectories(handles, side_str)
% ea_resultfig_update_trajectories(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, side_ids, ~] = ea_detsidestr(side_str);

resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
merhandles = getappdata(resultfig, 'merhandles');
mertoggles = getappdata(handles.mercontrolfig, 'mertoggles');

curr_fig_items = resultfig.CurrentAxes.Children;
tag_array = cell(1, length(curr_fig_items));
for cfi = 1:length(curr_fig_items)
    tag_array{cfi} = curr_fig_items(cfi).Tag;
end

set(0, 'CurrentFigure', resultfig); hold on;
for sid = side_ids
    for pos_ix = 1:length(merstruct.tract_info)
        pos_str = merstruct.tract_info(pos_ix).label;
        trajectory = merstruct.currentmer.(pos_str).trajectory{sid};
        h = merhandles.(pos_str){sid};
        if ~isempty(h)
            set(h, 'XData', trajectory(:,1)');
            set(h, 'YData', trajectory(:,2)');
            set(h, 'ZData', trajectory(:,3)');
        else
            merhandles.(pos_str){sid} = plot3(...
                trajectory(:,1), trajectory(:,2), trajectory(:,3),...
                'color', merstruct.tract_info(pos_ix).color, 'linew', 5,...
                'tag', [pos_str '_' side_strs{sid}]);
        end
        if mertoggles.togglestates(sid, pos_ix)
            set(merhandles.(pos_str){sid}, 'Visible', 'on');
        else
            set(merhandles.(pos_str){sid}, 'Visible', 'off');
        end
    end
end
setappdata(resultfig, 'merhandles', merhandles);