function ea_unified_perm_style_axes(ax)
% Applies a simple, publication-ready style to an axes: white background,
% no bounding box, no grid, outward ticks, clean sans-serif labels. Used
% consistently across all permutation-testing plots.

if ~exist('ax', 'var') || isempty(ax)
    ax = gca;
end

set(ancestor(ax, 'figure'), 'Color', 'w');
set(ax, 'Color', 'w');
box(ax, 'off');
grid(ax, 'off');
set(ax, 'TickDir', 'out');
set(ax, 'FontName', 'Helvetica', 'FontSize', 11);
set(ax, 'LineWidth', 0.75);
ax.XAxis.Color = [0.15, 0.15, 0.15];
ax.YAxis.Color = [0.15, 0.15, 0.15];
