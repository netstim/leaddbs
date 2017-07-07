function ea_resultfig_updatemarkers(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, ~, ~] = ea_detsidestr(side_str);

resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
mermarkers = getappdata(resultfig, 'mermarkers');

[shape.x, shape.y, shape.z] = sphere(20);  %20-voxel 3D sphere

set(0, 'CurrentFigure', resultfig); hold on;

% Go through each marker and plot it if not plotted already.
for marker_ix = 1:length(mermarkers)
    marker = mermarkers(marker_ix);
    
    if any(strcmpi(marker.side, side_strs))  && isempty(marker.handle)
        switch marker.markertype
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
        
        % Create text item
        marker_tag.depth = num2str(marker.depth);
        textpos = marker.coords_mm;
        textpos(1) = 3.2 * textpos(1) / abs(textpos(1)) + textpos(1);
        marker_tag.string = sprintf('%s.%s.%s: %smm',...
            marker.side(1), marker.tract(1:3), htag, marker_tag.depth);
        marker_tag.handle = text(textpos(1), textpos(2), textpos(3),...
            marker_tag.string, 'Color', 'w',...
            'HorizontalAlignment', 'center',...
            'Visible', merstruct.tag.visible);
        marker.tag = marker_tag;
        
        % Create sphere item
        marker_shape.x = shape.x*merstruct.markersize + marker.coords_mm(1);
        marker_shape.y = shape.y*merstruct.markersize + marker.coords_mm(2);
        marker_shape.z = shape.z*merstruct.markersize + marker.coords_mm(3);
        marker.handle = surf(marker_shape.x, marker_shape.y, marker_shape.z,...
            'FaceColor', fcolor, 'EdgeColor', 'none',...
            'FaceAlpha', 0.7, 'tag', htag);
        
        mermarkers(marker_ix) = marker;
    elseif isfield(marker.tag, 'handle') && ~isempty(marker.tag.handle)
        % Handle already exists. Just update its visibility.
        set(marker.tag.handle, 'Visible', merstruct.tag.visible);
    end
end
setappdata(resultfig, 'mermarkers', mermarkers);