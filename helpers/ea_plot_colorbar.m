function h = ea_plot_colorbar(cmap, width, orientation, titletxt, tick, ticklabel, target)
% Plot a standalone colorbar
%
% Parameters:
%     cmap: colormap to be plotted
%     width: width/height of the colorbar (vertical/horizontal mode)
%     orientation: 'v' for vertical or 'h' for horizontal
%     titletxt: tile of the colorbar figure
%     tick: tick locations
%     ticklabel: tick labels
%     target: figure or axes to plot the colorbar
%
% Output:
%     handle to the colorbar image figure
%
% Examples:
%     h = ea_plot_colorbar(colormap, 10, 'h', 'Cool');
%     ea_plot_colorbar(jet(128), [], 'v', '');

if ~exist('target', 'var')
    h = figure('Name', 'Colorbar');
    h.Position(3:4) = [370, 440];
    target = gca;
    plotInSeparateFigure = 1;
else
    plotInSeparateFigure = 0;
end

map = colormap(target, cmap);

if ~exist('width', 'var') || isempty(width)
    width = ceil(length(cmap)/16);
end

if ~exist('orientation', 'var') || isempty(orientation)
    orientation = 'v';
end

if ~exist('titletxt', 'var') || isempty(titletxt)
    titletxt = '';
end

switch lower(orientation)
    case {'v', 'vert', 'vertical'}
        image(target, repmat(cat(3, map(:,1), map(:,2), map(:,3)), 1, width));

        % Remove xticks
        set(gca, 'xtick', []);

        % Set yticks
        set(gca,'YAxisLocation','right')
        if exist('tick', 'var')
            set(gca, 'ytick', tick);
            if isempty(tick)
                set(target, 'yticklabel', []);
            elseif exist('ticklabel', 'var')
                if length(tick) == length(ticklabel)
                    set(target, 'yticklabel', ticklabel);
                else
                    error('tick and ticklabel should have the same length!');
                end
            else
                error('Please also specify ticklabel!');
            end
        else
            tick = get(target, 'ytick');
            tick = tick(mod(tick,1)==0);
            set(target, 'ytick', tick);
            set(target, 'yticklabel', tick);
        end

    case {'h', 'horz', 'horizontal'}
        image(target, repmat(cat(3, map(:,1)', map(:,2)', map(:,3)'), width, 1));

        % Remove yticks
        set(target, 'ytick', []);

        % Set xticks
        if exist('tick', 'var')
            set(target, 'xtick', tick);
            if exist('ticklabel', 'var')
                if length(tick) == length(ticklabel)
                    set(target, 'xticklabel', ticklabel);
                else
                    error('tick and ticklabel should have the same length!');
                end
            else
                error('Please also specify ticklabel!');
            end
        else
            tick = get(target, 'xtick');
            set(target, 'xtick', [0.5, tick]);
            set(target, 'xticklabel', [0, tick]);
        end

    otherwise
        error('Unknown colorbar orientation!');
end

set(target, 'Color', 'none');
set(get(target,'XLabel'), 'Visible', 'off');
set(get(target,'YLabel'), 'Visible', 'off');
axis(target, 'equal', 'tight', 'xy');

if plotInSeparateFigure
    title(titletxt);
end
