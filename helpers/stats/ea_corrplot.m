function [h,R,p,g]=ea_corrplot(X,Y,labels,corrtype,group1,group2,pperm,colors,markers)
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
% ea_corrplot(X,Y,{'Example Correlation','Age','Disease Duration'},'spearman',group1,group2)

if ~exist('labels','var')
    labels={'','X','Y'};
end

if ~(length(labels)==3) % assume only title provided
    labels{2}='X'; labels{3}='Y';
end

if ~exist('corrtype','var')
    corrtype='Pearson';
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

switch corrtype
    case {'permutation_spearman','permutation'}
        for tries=1:3
            try
                [R,p]=ea_permcorr(X,Y,'spearman');
            end
            if exist('R','var')
                break
            end
        end
    case 'permutation_pearson'
        for tries=1:3
            try
                [R,p]=ea_permcorr(X,Y,'pearson');
            end

            if exist('R','var')
                break
            end
        end
    otherwise
        [R,p]=corr(X,Y,'rows','pairwise','type',corrtype);
end

%Corner histogram
g=gramm('x',X,'y',Y);
if isempty(group1) && isempty(group2)
    g.geom_point();
    g.set_color_options(colorOptions{:});
else
    g.set_color_options('chroma',0,'lightness',30);
end
g.set_point_options('markers', markers, 'base_size', 7);
g.stat_glm();

pv=p;
pstr='p';
if exist('pperm','var') && ~isempty(pperm)
    pv=pperm;
    pstr='p(perm)';
end

if pv >= 0.001 % Show p = 0.XXX when p > 0.001
    pstr = [pstr, ' = ', sprintf('%.3f',pv)];
else
    % pstr = [pstr, ' = ', sprintf('%.1e',pv)]; % Show p = X.Xe-X
    signCheck=zeros(1,17);
    for i=1:length(signCheck)
        signCheck(i)=eval(['pv<1e-',num2str(i),';']);
    end
    pstr = [pstr, ' < 1e-', num2str(find(diff(signCheck),1))]; % Show p < 1e-X
end

g.set_title([labels{1}, ' [R = ', sprintf('%.2f',R), '; ', pstr, ']'], 'FontSize', 20);
g.set_names('x',labels{2},'y',labels{3});
g.set_text_options('base_size',22);
g.no_legend();

ratio = 7/8;
Width = 550;
Height = Width*ratio;
h=figure('Position',[100 100 Width Height]);
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
