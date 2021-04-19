function [h,pv,Rsquared,F,mdl,bb,dev,allstats]=ea_glmplot(X,y,labels,distribution,group1,group2,colors,markers,switchxy)

if ~(size(y,2)==1)
    ea_warning('Assuming X and y were switched. Switching variables.');
    Xn=y;
    y=X;
    X=Xn;
end

if ~exist('labels','var')
    labels={'','X','Y'};
end

if ~(length(labels)==3) % assume only title provided
    labels{2}='X'; labels{3}='Y';
end

if ~exist('corrtype','var')
    distribution='normal';
end

if ~exist('group1','var')
    group1=[];
else
    if ~isstruct(group1)
        group1s.idx=group1;
        group1s.tag='color';
        clear group1
        group1=group1s;
    end
end

if ~exist('group2','var')
    group2=[];
else
    if ~isstruct(group2)
        group2s.idx=group2;
        group2s.tag='color';
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

if exist('markers', 'var') && ~isempty(markers)
    if ~isempty(group2)
        if numel(markers) < numel(unique(group2.idx))
            error('Number of custom markers is less than number of categories!');
        end
    end
else
    markers = {'o' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '*' '+' 'x'};
end

if ischar(map)
    colorOptions = {'map', map};
else
    colorOptions = {'map', map, 'n_color', size(map,1), 'n_lightness', 1};
end

mdl=fitglm(X,y,'distribution',distribution);
[bb,dev,allstats]=glmfit(X,y,distribution);
yhat=predict(mdl,X);

if exist('switchxy','var') && switchxy
    g=gramm('x',y,'y',yhat); % data needs to be put in "reversed" for gramm.
else
    g=gramm('x',yhat,'y',y); % data needs to be put in "reversed" for gramm.
end

if isempty(group1) && isempty(group2)
    g.geom_point();
    g.set_color_options(colorOptions{:});
else
    g.set_color_options('chroma',0,'lightness',30);
end

g.set_point_options('markers', markers, 'base_size', 7);
g.stat_glm('distribution',distribution,'fullrange','false','fullrange','false');

[~,~,~,~,stats]=regress(y,ea_addone(X));
%g.geom_abline();
%g.stat_cornerhist('edges',[ea_nanmean(X)-1*ea_nanstd(X):0.2:ea_nanmean(X)+1*ea_nanstd(X)],'aspect',0.6,'location',max(X));
pstr='p';
pv=stats(3);
F=stats(2);
Rsquared=mdl.Rsquared.Ordinary;

g.set_title([labels{1},' [R2 = ',sprintf('%.2f',Rsquared),'; F-stat = ',sprintf('%.2f',F),'; ',pstr,' = ',sprintf('%.3f',pv),']'],'FontSize',20);
g.set_names('x',labels{2},'y',labels{3});
g.set_text_options('base_size',22);
g.no_legend();

h=figure('Position',[100 100 550 550]);
g.draw();

if ~isempty(group2) && ~isempty(group1)
    g.update('marker',group2.idx,'color',group1.idx);
    g.set_color_options(colorOptions{:});
    g.set_names('marker',group2.tag,'color',group1.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
    set(h,'Position',[100 100 650 550]);
elseif ~isempty(group2) && isempty(group1)
    g.update('marker',group2.idx);
    g.set_color_options();
    g.set_names('marker',group2.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
    set(h,'Position',[100 100 650 550]);
elseif isempty(group2) && ~isempty(group1)
    g.update('color',group1.idx);
    g.set_color_options(colorOptions{:});
    g.set_names('color',group1.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
    set(h,'Position',[100 100 650 550]);
end

set([g.results.geom_point_handle],'MarkerSize',7);
set([g.results.geom_point_handle],'MarkerEdgeColor','w');
