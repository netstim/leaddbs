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

% Synchronize markers and handles.

% Get a key for each handle and for each marker in merstruct
uq_sides = unique(cat(2, {merstruct.Markers.side}, {merhandles.markers.side}));
uq_labels = unique(cat(2, {merstruct.Markers.tract_label}, {merhandles.markers.label}));
[~, h_side_keys] = ismember({merhandles.markers.side}, uq_sides);
[~, h_label_keys] = ismember({merhandles.markers.label}, uq_labels);
h_keys = [h_side_keys' h_label_keys' [merhandles.markers.depth]'];
[~, s_side_keys] = ismember({merstruct.Markers.side}, uq_sides);
[~, s_label_keys] = ismember({merstruct.Markers.tract_label}, uq_labels);
s_keys = [s_side_keys' s_label_keys' [merstruct.Markers.depth]'];

% Clear handles without matching keys in merstruct. e.g., after undo
if ~isempty(h_keys)
    [~, ia] = setdiff(h_keys, s_keys, 'rows');
    for h_ix = 1:length(ia)
        h_id = ia(h_ix);
        delete(merhandles.markers(h_id).h);
        delete(merhandles.markers(h_id).tag);
    end
    merhandles.markers(ia) = [];
    
    % Update visibility of handles that do have matching keys in merstruct
    ia = find(ismember(h_keys, s_keys, 'rows'));
    for h_ix = 1:length(ia)
        h_id = ia(h_ix);
        hm = merhandles.markers(h_id);
        if isfield(hm, 'h') && isvalid(hm.h)
            set(hm.tag, 'Visible', merstruct.Config.vis.tag_visible);
            set(hm.h, 'Visible', 'on');
        end
    end
    
    % Create new handles for merstruct markers without keys in merhandles.
    [~, ia] = setdiff(s_keys, h_keys, 'rows');
else
    ia = 1:size(s_keys, 1);
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
    textpos = marker.coords.(space_type);
    textpos(1) = 3.2 * textpos(1) / abs(textpos(1)) + textpos(1);
    tag_string = sprintf('%s.%s.%s: %.2fmm',...
        marker.side(1), marker.tract_label(1:3), htag(1:3), marker.depth);
    tag_handle = text(textpos(1), textpos(2), textpos(3),...
        tag_string, 'Color', 'w',...
        'HorizontalAlignment', 'center',...
        'Visible', merstruct.Config.vis.tag_visible);
    new_hmark.tag = tag_handle;

    % Create sphere item
    marker_shape.x = shape.x*merstruct.Config.vis.markersize + marker.coords.(space_type)(1);
    marker_shape.y = shape.y*merstruct.Config.vis.markersize + marker.coords.(space_type)(2);
    marker_shape.z = shape.z*merstruct.Config.vis.markersize + marker.coords.(space_type)(3);
    new_hmark.h = surf(marker_shape.x, marker_shape.y, marker_shape.z,...
        'FaceColor', fcolor, 'EdgeColor', 'none',...
        'FaceAlpha', 0.7, 'tag', htag);

    merhandles.markers(m_id) = new_hmark;
end
setappdata(handles.mercontrolfig, 'merhandles', merhandles);