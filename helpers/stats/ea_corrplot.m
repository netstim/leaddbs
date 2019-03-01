function [h,R,p]=ea_corrplot(X,Y,labels,corrtype,group1,group2,pperm)
% Wrapper for gramm to produce a simple correlation plot. Group1 denotes
% colors, Group2 Markers.
% (c) Andreas Horn 2019 Charite Berlin

% Example usage:
% ---------------
% X=randn(100,1);
% Y=X.*randn(100,1);
% 
% ea_corrplot(X,Y)
% 
% group1cell={'Prague','Berlin','London','Moscow','Paris','Madrid'};
% group1.idx=group1cell(ceil(rand(100,1)*5));
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
end


if ~exist('group2','var')
    group2=[];
end

switch corrtype
    case {'permutation_spearman','permutation'}
        [R,p]=ea_permcorr(X,Y,'spearman');
    case 'permutation_pearson'
        [R,p]=ea_permcorr(X,Y,'pearson');
    otherwise
        [R,p]=corr(X,Y,'rows','pairwise','type',corrtype);
end


%Corner histogram
g=gramm('x',X,'y',Y);
if isempty(group1) && isempty(group2)
    g.geom_point();
else
    g.set_color_options('chroma',0,'lightness',30);
end
g.stat_glm();
%g.geom_abline();
%g.stat_cornerhist('edges',[ea_nanmean(X)-1*ea_nanstd(X):0.2:ea_nanmean(X)+1*ea_nanstd(X)],'aspect',0.6,'location',max(X));
pstr='p';
pv=p;
if exist('pperm','var')
    pv=pperm;
    pstr='p(perm)';
end

g.set_title([labels{1},' [R = ',sprintf('%.2f',R),'; ',pstr,' = ',sprintf('%.3f',pv),']'],'FontSize',20);
g.set_names('x',labels{2},'y',labels{3});
g.set_text_options('base_size',22);
g.no_legend();
h=figure('Position',[100 100 550 550]);
g.draw();
if ~isempty(group2) && ~isempty(group1)
    g.update('marker',group2.idx,'color',group1.idx);
    g.set_color_options();
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
elseif ~isempty(group1) && ~isempty(group2)
    g.update('color',group1.idx);
    g.set_color_options();
    g.set_names('color',group1.tag,'x',labels{2},'y',labels{3});
    g.geom_point();
    g.draw();
    set(h,'Position',[100 100 650 550]);
end

set([g.results.geom_point_handle],'MarkerSize',7);
set([g.results.geom_point_handle],'MarkerEdgeColor','w');

