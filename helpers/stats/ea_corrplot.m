function [h,R,p,g]=ea_corrplot(X,Y,permutation,labels,group1,group2,colors,markers)
% Wrapper for gramm to produce a simple correlation plot.
% Group1 denotes colors, Group2 Markers.
% Can also specify custom colors
% (c) Andreas Horn 2019 Charite Berlin

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
% ea_corrplot(X,Y,1000,{'Example Correlation','Age','Disease Duration'},group1,group2)

if isrow(X)
    X = X';
end

if isrow(Y)
    Y = Y';
end

if ~exist('permutation', 'var')
    permutation = 0;
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
    labels{4} = 'Default corrplot';
end

if ~exist('group1','var')
    group1=[];
else
    if ~isempty(group1) && ~isstruct(group1)
        group1s.idx=group1;
        group1s.tag='group';
        clear group1
        group1=group1s;
    end
end

if ~exist('group2','var')
    group2=[];
else
    if ~isempty(group2) && ~isstruct(group2)
        group2s.idx=group2;
        group2s.tag='type';
        clear group2
        group2=group2s;
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

if contains(labels{4}, 'nested LOO', 'IgnoreCase', true)
    g = gramm('x', X, 'y', Y, 'color', group1.idx);
else
    g = gramm('x', X, 'y', Y);
    if isempty(group1) && isempty(group2)
        g.geom_point();
        g.set_color_options(colorOptions{:});
    else
        g.set_color_options('chroma', 0, 'lightness', 30);
    end
end

g.set_point_options('markers', markers, 'base_size', 7);
g.stat_glm();

if isnumeric(permutation)
    if ~permutation
        [R_linear, p_linear] = corr(X, Y, 'rows', 'pairwise', 'type', 'Pearson');
        [R_rank, p_rank] = corr(X, Y, 'rows', 'pairwise', 'type', 'Spearman');
        pstr_linear = getPstr(p_linear, 'p');
        pstr_rank = getPstr(p_rank, 'p');
    else
        [R_linear, p_linear] = ea_permcorr(X, Y, 'Pearson', permutation);
        [R_rank, p_rank] = ea_permcorr(X, Y, 'Spearman', permutation);
        pstr_linear = getPstr(p_linear, 'p (perm)');
        pstr_rank = getPstr(p_rank, 'p (perm)');
    end
elseif ischar(permutation)
    switch permutation
        case {'no', 'noperm', 'nopermutation'}
            [R_linear, p_linear] = corr(X, Y, 'rows', 'pairwise', 'type', 'Pearson');
            [R_rank, p_rank] = corr(X, Y, 'rows', 'pairwise', 'type', 'Spearman');
            pstr_linear = getPstr(p_linear, 'p');
            pstr_rank = getPstr(p_rank, 'p');
        case {'yes', 'perm', 'permutation'}
            [R_linear, p_linear] = ea_permcorr(X, Y, 'Pearson');
            [R_rank, p_rank] = ea_permcorr(X, Y, 'Spearman');
            pstr_linear = getPstr(p_linear, 'p (perm)');
            pstr_rank = getPstr(p_rank, 'p (perm)');
    end
end

if contains(labels{4}, 'nested LOO', 'IgnoreCase', true)
    g.set_title({'Mean and STD of linear models from nested LOO', ['Slope: ',labels{5}], ['Intercept: ',labels{6}]})
elseif length(labels) == 4
    g.set_title({[labels{1}], ['Spearman: [R = ', sprintf('%.2f',R_rank), '; ', pstr_rank, ']'], ['Pearson: [R = ', sprintf('%.2f',R_linear), '; ', pstr_linear, ']']});
elseif length(labels) > 4
    g.set_title({[labels{1}], ['Spearman: [R = ', sprintf('%.2f',R_rank), '; ', pstr_rank, ']'], ['Pearson: [R = ', sprintf('%.2f',R_linear), '; ', pstr_linear, ']'], [labels{5}, '; ',labels{6}, '; ',labels{7}]});
end


fs=3 * (25/length(labels{1}));
if fs>30
    fs=40;
end
% g.set_title([labels{1}, ' [R = ', sprintf('%.2f',R), '; ', pstr, ']'], 'FontSize', fs);

g.set_names('x',labels{2},'y',labels{3});
fs=2*(35/max(cellfun(@length,labels(2:3))));
if fs>30
    fs=40;
end
% g.set_text_options('base_size',fs);

g.no_legend();

ratio = 7/8;
Width = 550;
Height = Width*ratio;
h=figure('Name',labels{4},'NumberTitle','off','Position',[100 100 Width Height]);
g.draw();

gtitle = g.title_axe_handle.Children;
gtitle.Units = 'pixels';

if gtitle.Extent(3) > Width
    % Shift the title a bit
    if isempty(group2) && isempty(group1)
        gtitle.Position(1) = 376;
    else
        gtitle.Position(1) = 246;
    end

    % Calculate new figure size
    Width = gtitle.Extent(3)+100;
    Height = Width*ratio;
end

if ~isempty(group2) && ~isempty(group1)
    g.update('marker',group2.idx,'color',group1.idx);
    g.set_color_options(colorOptions{:});
    g.set_names('marker',group2.tag,'color',group1.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
elseif ~isempty(group2) && isempty(group1)
    g.update('marker',group2.idx);
    g.set_color_options();
    g.set_names('marker',group2.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
elseif isempty(group2) && ~isempty(group1)
    g.update('color',group1.idx);
    g.set_color_options(colorOptions{:});
    g.set_names('color',group1.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
end

set(h,'Position',[100 100 Width Height]);
set([g.results.geom_point_handle],'MarkerSize',7);
set([g.results.geom_point_handle],'MarkerEdgeColor','w');


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
