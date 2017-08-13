function ea_resultfig_updatemarkers(handles)

resultfig = getappdata(handles.mercontrolfig, 'resultfig');
options = getappdata(handles.mercontrolfig, 'options');
merstruct = getappdata(handles.mercontrolfig, 'merstruct');
merhandles = getappdata(handles.mercontrolfig, 'merhandles');

if options.native
    space_type = 'native';
else
    space_type = 'mni';
end

[shape.x, shape.y, shape.z] = sphere(20);  %20-voxel 3D sphere

set(0, 'CurrentFigure', resultfig); hold on;

% Clear handles without matching keys in merstruct. e.g., after undo
[s_keys, h_keys] = get_keys_for_markers_handles(merstruct.Markers, merhandles.markers);
if ~isempty(h_keys)
    if isempty(s_keys)
        ia = intersect(1:size(merhandles.markers,2),1:length(h_keys));
    else
        [~, ia] = setdiff(h_keys, s_keys, 'rows');
    end
    for h_ix = 1:length(ia)
        h_id = ia(h_ix);
        try delete(merhandles.markers(h_id).h); end
        try delete(merhandles.markers(h_id).tag); end
    end
    merhandles.markers(ia) = [];
end
    
% Update visibility of handles that do have matching keys in merstruct
[s_keys, h_keys] = get_keys_for_markers_handles(merstruct.Markers, merhandles.markers);
if ~isempty(s_keys) && ~isempty(h_keys)
    ia = find(ismember(h_keys, s_keys, 'rows'));
    for h_ix = 1:length(ia)
        h_id = ia(h_ix);
        hm = merhandles.markers(h_id);
        if isfield(hm, 'h') && isvalid(hm.h)
            set(hm.tag, 'Visible', merstruct.Config.vis.tag_visible);
            set(hm.h, 'Visible', 'on');
        end
    end
end
    
% Create new handles for merstruct markers without keys in merhandles.
if isempty(h_keys)
    ia = 1:size(s_keys, 1);  % All markers
else
    [~, ia] = setdiff(s_keys, h_keys, 'rows');  % Missing markers only.
end
for m_ix = 1:length(ia)
    m_id = ia(m_ix);
    marker = merstruct.Markers(m_id);
    switch marker.type
        case 'Generic'
            htag = 'Generic';
            fcolor = [0.5 0.5 0];
        case 'MER recording'
            htag = 'MER';
            fcolor = [0.5 0 0];
        case 'LFP recording'
            htag = 'LFP';
            fcolor = [0 0.5 0];
        case 'Top border'
            htag = 'Top';
            fcolor = [0 0 0.5];
        case 'Bottom border'
            htag = 'Bottom';
            fcolor = [0 0 0.5];
    end
    new_hmark = struct('side', marker.side, 'label', marker.tract_label,...
        'depth', marker.depth, 'h', [], 'tag', []);
    % Create text item
    marker_coords = merstruct.getMarkerPosition(marker, space_type);
    textpos = marker_coords;
    textpos(1) = 3.2 * textpos(1) / abs(textpos(1)) + textpos(1);
    tag_string = sprintf('%s.%s.%s: %.2fmm',...
        marker.side(1), marker.tract_label(1:3), htag(1:3), marker.depth);
    tag_handle = text(textpos(1), textpos(2), textpos(3),...
        tag_string, 'Color', 'w',...
        'HorizontalAlignment', 'center',...
        'Visible', merstruct.Config.vis.tag_visible);
    new_hmark.tag = tag_handle;

    % Create sphere item
    marker_shape.x = shape.x*merstruct.Config.vis.markersize + marker_coords(1);
    marker_shape.y = shape.y*merstruct.Config.vis.markersize + marker_coords(2);
    marker_shape.z = shape.z*merstruct.Config.vis.markersize + marker_coords(3);
    new_hmark.h = surf(marker_shape.x, marker_shape.y, marker_shape.z,...
        'FaceColor', fcolor, 'EdgeColor', 'none',...
        'FaceAlpha', 0.7, 'tag', htag);

    merhandles.markers(m_id) = new_hmark;
end
setappdata(handles.mercontrolfig, 'merhandles', merhandles);

function [s_keys, h_keys] = get_keys_for_markers_handles(s_markers, h_markers)
% Get a key for each marker in merstruct and for each marker handle.
uq_sides = unique(cat(2, {s_markers.side}, {h_markers.side}));
uq_labels = unique(cat(2, {s_markers.tract_label}, {h_markers.label}));
[~, h_side_keys] = ismember({h_markers.side}, uq_sides);
[~, h_label_keys] = ismember({h_markers.label}, uq_labels);
h_keys = [h_side_keys' h_label_keys' [h_markers.depth]'];
[~, s_side_keys] = ismember({s_markers.side}, uq_sides);
[~, s_label_keys] = ismember({s_markers.tract_label}, uq_labels);
s_keys = [s_side_keys' s_label_keys' [s_markers.depth]'];