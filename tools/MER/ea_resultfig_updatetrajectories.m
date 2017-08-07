function ea_resultfig_updatetrajectories(handles, side_str)
% ea_resultfig_update_trajectories(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, ~, ~] = ea_detsidestr(side_str);

resultfig = getappdata(handles.mercontrolfig, 'resultfig');
options = getappdata(handles.mercontrolfig, 'options');
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merhandles = getappdata(handles.mercontrolfig, 'merhandles');

spc = 'mni';
if options.native
    spc = 'native';
end

uqlabels = unique({merstruct.MERTrajectories.label}, 'stable');

set(0, 'CurrentFigure', resultfig); hold on;
for traj_ix = 1:length(merstruct.MERTrajectories)
    mertraj = merstruct.MERTrajectories(traj_ix);
    if any(strcmpi(side_strs, mertraj.side))
        bH = strcmpi({merhandles.traj.side}, mertraj.side) ...
            & strcmpi({merhandles.traj.label}, mertraj.label);
        h = [merhandles.traj(bH).h];
        coords = merstruct.getMERTrajectory(mertraj, spc);
        if length(h) == 1
            set(h, 'XData', coords(:, 1)');
            set(h, 'YData', coords(:, 2)');
            set(h, 'ZData', coords(:, 3)');
        else
            delete(h);
            merhandles.traj(bH) = [];
            
            this_color = merhandles.color_list(strcmpi(uqlabels, mertraj.label), :);
            merhandles.traj(end + 1) = struct('side', mertraj.side,...
                'label', mertraj.label,...
                'h', plot3(coords(:, 1), coords(:, 2), coords(:, 3),...
                'color', this_color, 'linew', 5,...
                'tag', [mertraj.label '_' mertraj.side]));
            bH = false(size(merhandles.traj)); bH(end) = true;
        end
        bToggle = strcmpi({merstruct.Toggles.togglestates.side}, mertraj.side) ...
            & strcmpi({merstruct.Toggles.togglestates.label}, mertraj.label);
        if merstruct.Toggles.togglestates(bToggle).value
            set(merhandles.traj(bH).h, 'Visible', 'on');
        else
            set(merhandles.traj(bH).h, 'Visible', 'off');
        end
    end
end
setappdata(handles.mercontrolfig, 'merhandles', merhandles);