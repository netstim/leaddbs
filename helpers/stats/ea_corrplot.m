function [h,R,p,g] = ea_corrplot(X,Y,permutation,labels,group1,group2,colors,markers,plottype,h)
% Wrapper for gramm to create a correlation plot while showing stats based
% on [permuted] rank and linear corrlation in the title.
%
%   permutation can be:
%       numeric: integer (larger than 1) for number of permutation
%                float for external provided p value (for example, p
%                value can be calculated externally based on permuted R-Map
%                or fiber tracts).
%       char:    'no', 'noperm', or 'nopermutation' for no permutation.
%                'yes', 'perm', or 'permutation' for permutation (default
%                number of permutation is 5000 as defined in ea_permcorr).
%       logical: false for no permutation
%                true for permutation (default number of permutation is
%                5000 as defined in ea_permcorr)
%
%   labels: labels{1} for figure title
%           labels{2} for xlabel
%           labels{3} for ylabel
%           labels{4} for figure name
%           labels{5:end} for other info appended to the figure title
%
%   group1: color group
%
%   group2: marker group
%
%   colors: custom colors
%
%   markers: custom markers
%
%   plottype: 'linear' (default) or 'rank' plot
%
%   (c) Andreas Horn 2019 Charite Berlin

% Example usage:
% ---------------
% X=randn(100,1);
% Y=X.*randn(100,1);
%
% ea_corrplot(X,Y)
%
% group1cell={'Prague','Berlin','London','Moscow','Paris','Madrid'};
% group1.idx=group1cell(ceil(rand(100,1)*6));
% group1.tag='Cohort';
%
% group2cell={'Parkinson','Alzheimer'};
% group2.idx=group2cell(ceil(rand(100,1)*2));
% group2.tag='Disease';
%
% ea_corrplot(X,Y,5000,{'Example Correlation','Age','Disease Duration'},group1,group2)

if isrow(X)
    X = X';
end

if isrow(Y)
    Y = Y';
end

if ~exist('permutation', 'var')
    permutation = nan;
end

if ~exist('labels', 'var')
    labels={'', 'X', 'Y'};
end

if length(labels) < 2
    labels{2} = 'X';
end

if length(labels) < 3
    labels{3} = 'Y';
end

if length(labels) < 4
    labels{4} = 'corrplot';
end

if ~exist('group1','var')
    group1 = [];
else
    if ~isempty(group1) && ~isstruct(group1)
        group1s.idx = group1;
        group1s.tag = 'group';
        clear group1
        group1 = group1s;
    end
end

if ~exist('group2','var')
    group2 = [];
else
    if ~isempty(group2) && ~isstruct(group2)
        group2s.idx = group2;
        group2s.tag = 'type';
        clear group2
        group2 = group2s;
    end
end

if exist('colors', 'var') && ~isempty(colors)
    map = colors;
    if isnumeric(map) && ~isempty(group1)
        if size(map,1) ~= numel(unique(group1.idx))
            error('Number of custom colors doesn''t match number of categories!');
        end
    end
else
    map = 'lch';
end

if ischar(map)
    colorOptions = {'map', map};
else
    colorOptions = {'map', map, 'n_color', size(map,1), 'n_lightness', 1};
end

if exist('markers', 'var') && ~isempty(markers)
    if ~isempty(group2)
        if numel(markers) < numel(unique(group2.idx))
            error('Number of custom markers is less than number of categories!');
        end
    end
else
    markers = {'o' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '*' '+' 'x'};
end

if ~exist('plottype', 'var')
    plottype = 'linear';
else
    switch lower(plottype)
        case {'linear', 'pearson'}
            plottype = 'linear';
        case {'rank', 'spearman'}
            plottype = 'rank';
    end
end

ratio = 7/8;
Width = 550;
Height = Width*ratio;
if ~exist('h', 'var')
    h = figure('Name', [labels{4}, ' (', plottype, ')'] , 'NumberTitle', 'off', 'Position', [100 100 Width Height]);
    inputData = {X, Y, permutation, labels, group1, group2, map, markers};
    setappdata(h, 'inputData', inputData);
    setappdata(h, 'plottype', plottype);
    toolbar = findall(h, 'Tag', 'FigureToolBar');
    uipushtool(toolbar, 'Tag', 'TogglePlot', 'CData', ea_get_icn('corrplot'), 'ClickedCallback', {@(src, evt) toggleplot(h)});
end

toggle = findall(h, 'Tag', 'TogglePlot');
plottype = getappdata(h, 'plottype');
if strcmp(plottype, 'rank')
    Xplot = ea_vec2rank(X);
    Yplot = ea_vec2rank(Y);
    h.Name = strrep(h.Name, '(linear)', '(rank)');
    toggle.Tooltip = 'Switch to Linear Plot';
else
    Xplot = X;
    Yplot = Y;
    h.Name = strrep(h.Name, '(rank)', '(linear)');
    toggle.Tooltip = 'Switch to Rank Plot';
end

if contains(labels{4}, 'nested LOO', 'IgnoreCase', true)
    g = gramm('x', Xplot, 'y', Yplot, 'color', group1.idx);
else
    g = gramm('x', Xplot, 'y', Yplot);
    if isempty(group1) && isempty(group2)
        g.geom_point();
        g.set_color_options(colorOptions{:});
    else
        g.set_color_options('chroma', 0, 'lightness', 30);
    end
end

g.set_point_options('markers', markers, 'base_size', 7);

g.stat_glm();

g.set_names('x', labels{2}, 'y', labels{3});

baseFontSize = 12;
g.set_text_options('base_size', baseFontSize);

R_linear = getappdata(h, 'R_linear');
p_linear = getappdata(h, 'p_linear');
pstr_linear = getappdata(h, 'pstr_linear');
R_rank = getappdata(h, 'R_rank');
p_rank = getappdata(h, 'p_rank');
pstr_rank = getappdata(h, 'pstr_rank');

if isempty(R_linear)
    if isnumeric(permutation)
        if permutation > 1 && mod(permutation, 1) == 0 % integer, number of permutation
            [R_linear, p_linear] = ea_permcorr(X, Y, 'Pearson', permutation);
            [R_rank, p_rank] = ea_permcorr(X, Y, 'Spearman', permutation);
            pstr_linear = getPstr(p_linear, 'p (perm)');
            pstr_rank = getPstr(p_rank, 'p (perm)');
        else % external p-value provided (e.g., based on permuted R-Map or fiber tracts)
            [R_linear, p_linear] = ea_permcorr(X, Y, 'Pearson');
            [R_rank, p_rank] = ea_permcorr(X, Y, 'Spearman');
            pstr_linear = getPstr(p_linear, 'p (perm)');
            pstr_rank = getPstr(p_rank, 'p (perm)');
        end
    elseif ischar(permutation) || islogical(permutation)
        switch permutation % no permutation
            case {'no', 'noperm', 'nopermutation', false}
                [R_linear, p_linear] = corr(X, Y, 'rows', 'pairwise', 'type', 'Pearson');
                [R_rank, p_rank] = corr(X, Y, 'rows', 'pairwise', 'type', 'Spearman');
                pstr_linear = getPstr(p_linear, 'p');
                pstr_rank = getPstr(p_rank, 'p');
            case {'yes', 'perm', 'permutation', true} % default number of permutation defined in ea_permcorr
                [R_linear, p_linear] = ea_permcorr(X, Y, 'Pearson');
                [R_rank, p_rank] = ea_permcorr(X, Y, 'Spearman');
                pstr_linear = getPstr(p_linear, 'p (perm)');
                pstr_rank = getPstr(p_rank, 'p (perm)');
        end
    end
    setappdata(h, 'R_linear', R_linear);
    setappdata(h, 'p_linear', p_linear);
    setappdata(h, 'pstr_linear', pstr_linear);
    setappdata(h, 'R_rank', R_rank);
    setappdata(h, 'p_rank', p_rank);
    setappdata(h, 'pstr_rank', pstr_rank);
end

% external p-value provided (e.g., based on permuted R-Map or fiber tracts)
if isnumeric(permutation) && ~isnan(permutation) && permutation <= 1
    p_external = permutation;
    labels = [labels, getPstr(p_external, 'p (external)')];
end

R.spearman = R_rank;
R.pearson = R_linear;
p.spearman = p_rank;
p.pearson = p_linear;
if exist('p_external', 'var')
    p.external = p_external;
end

label_linear = ['Pearson: [R = ', sprintf('%.2f',R_linear), '; ', pstr_linear, ']'];
label_rank = ['Spearman: [R = ', sprintf('%.2f',R_rank), '; ', pstr_rank, ']'];

titleFontSize = 14;

if contains(labels{4}, 'nested LOO', 'IgnoreCase', true)
    g.set_title({'Mean and STD of linear models from nested LOO', ['Slope: ',labels{5}], ['Intercept: ',labels{6}]}, 'FontSize', titleFontSize)
elseif length(labels) == 4
    g.set_title({labels{1}, label_rank, label_linear}, 'FontSize', titleFontSize);
elseif length(labels) > 4
    g.set_title({labels{1}, label_rank, label_linear, strjoin(labels(5:end), '; ')}, 'FontSize', titleFontSize);
end

g.no_legend();

g.draw();

gtitle = g.title_axe_handle.Children;
gtitle.Units = 'pixels';

% Adapt figure size when title is oversized
if gtitle.Extent(3) > Width
    % Calculate new figure size
    Width = gtitle.Extent(3) + 100;
    Height = Width*ratio;
    % Make sure to set height first and then width, so the figure is
    % enlarged proportionally with the title position shifted accordingly.
    h.Position(4) = Height;
    h.Position(3) = Width;
end

if ~isempty(group2) && ~isempty(group1)
    g.update('marker', group2.idx, 'color', group1.idx);
    g.set_color_options(colorOptions{:});
    g.set_names('marker', group2.tag, 'color', group1.tag, 'x', labels{2}, 'y', labels{3});
    g.geom_point();
    g.draw();
elseif ~isempty(group2) && isempty(group1)
    g.update('marker', group2.idx);
    g.set_color_options();
    g.set_names('marker', group2.tag, 'x', labels{2}, 'y', labels{3});
    g.geom_point();
    g.draw();
elseif isempty(group2) && ~isempty(group1)
    g.update('color', group1.idx);
    g.set_color_options(colorOptions{:});
    g.set_names('color', group1.tag, 'x', labels{2}, 'y', labels{3});
    g.geom_point();
    g.draw();
end

set(g.results.geom_point_handle, 'MarkerSize', 7);
set(g.results.geom_point_handle, 'MarkerEdgeColor', 'w');


function pstr = getPstr(p, prefix)
if p >= 0.001 % Show p = 0.XXX when p >= 0.001
    pstr = [prefix, ' = ', sprintf('%.3f',p)];
else
    % pstr = [pstr, ' = ', sprintf('%.1e',pv)]; % Show p = X.Xe-X
    signCheck=zeros(1,16);
    for i=1:length(signCheck)
        signCheck(i) = eval(['p < 1e-',num2str(i),';']);
    end
    if all(signCheck)
        pstr = [prefix, ' < 1e-16']; % Show p < 1e-16
    else
        pstr = [prefix, ' < 1e-', num2str(find(diff(signCheck),1))]; % Show p < 1e-X
    end
end


function toggleplot(h)
clf(h);
inputData = getappdata(h, 'inputData');
plottype = getappdata(h, 'plottype');
if strcmp(plottype, 'linear')
    setappdata(h, 'plottype', 'rank');
elseif strcmp(plottype, 'rank')
    setappdata(h, 'plottype', 'linear');
end

ea_corrplot(inputData{:}, plottype, h);
