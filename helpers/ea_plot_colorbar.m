function h = ea_plot_colorbar(cmap, width, orientation, titletxt, tick, ticklabel)
% Plot a standalone colorbar
%
% Parameters:
%     cmap: colormap to be plotted
%     width: width/height of the colorbar (vertical/horizontal mode)
%     orientation: 'v' for vertical or 'h' for horizontal
%     title: tile of the colorbar
%     tick: tick locations
%     label: tick labels
%
% Output:
%     handle to the colorbar image figure
%
% Examples:
%     h = ea_plot_colorbar(colormap, 10, 'h', 'Cool');
%     ea_plot_colorbar(jet(128), [], 'v', '');

h = figure('Name', 'Colorbar');
h.Position(3:4) = [370, 440];

map = colormap(gcf, cmap);

if ~exist('width', 'var') || isempty(width)
    width = ceil(length(cmap)/16);
end

if ~exist('orientation', 'var') || isempty(orientation)
    orientation = 'v';
end

if ~exist('titletxt', 'var')
    titletxt = '';
end

switch lower(orientation)
    case {'v', 'vert', 'vertical'}
        image(repmat(cat(3, map(:,1), map(:,2), map(:,3)), 1, width));

        % Remove xticks
        set(gca, 'xtick', []);

        % Set yticks
        set(gca,'YAxisLocation','right')
        if exist('tick', 'var')
            set(gca, 'ytick', tick);
            if exist('ticklabel', 'var')
                if length(tick) == length(ticklabel)
                    set(gca, 'yticklabel', ticklabel);
                else
                    error('tick and ticklabel should have the same length!');
                end
            else
                error('Please also specify ticklabel!');
            end
        else
            tick = get(gca, 'ytick');
            set(gca, 'ytick', [0.5, tick]);
            set(gca, 'yticklabel', [0, tick]);
        end

    case {'h', 'horz', 'horizontal'}
        h = image(repmat(cat(3, map(:,1)', map(:,2)', map(:,3)'), width, 1));

        % Remove yticks
        set(gca, 'ytick', []);

        % Set xticks
        if exist('tick', 'var')
            set(gca, 'xtick', tick);
            if exist('ticklabel', 'var')
                if length(tick) == length(ticklabel)
                    set(gca, 'xticklabel', ticklabel);
                else
                    error('tick and ticklabel should have the same length!');
                end
            else
                error('Please also specify ticklabel!');
            end
        else
            tick = get(gca, 'xtick');
            set(gca, 'xtick', [0.5, tick]);
            set(gca, 'xticklabel', [0, tick]);
        end

    otherwise
        error('Unknown colorbar orientation!');
end

% Set up the axis
title(titletxt)
axis equal
axis tight
axis xy
