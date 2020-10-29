function [R_upd,p_upd,R,p,f]=ea_corrplot_col(X,Y,labels,corrtype,groupcolor,groupfiducial,pperm)
% this simple function is a small wrapper for a corrplot figure. Same
% interface as ea_corrplot

description=labels{1};
labels=labels(2:3);

% if nargin>4
%     if ~isempty(varargin{5})
%         if ischar(varargin{5}) || isequal([1 3],size(varargin{5}))
%             color=varargin{5};
%         else
%             groups=varargin{5};
%         end
%     else
%         groups=ones(length(X),1);
%     end
% else
%     groups=ones(length(X),1);
% end

if ~exist('corrtype','var')
corrtype='Pearson';
end

disp(description);

% support for nans
nnans=isnan(sum([X,Y],2));
X(nnans,:)=[];
Y(nnans,:)=[];

try groupcolor(nnans,:)=[]; end
try groupfiducial(nnans,:)=[]; end

switch corrtype
    case {'permutation_spearman','permutation'}
        [R,p]=ea_permcorr(X,Y,'spearman');
        R_upd=R; p_upd=p;
    case 'permutation_pearson'
        [R,p]=ea_permcorr(X,Y,'pearson');
        R_upd=R; p_upd=p;
    otherwise
        [R,p]=corr(X,Y,'rows','pairwise','type',corrtype);
        R_upd=R(2:end,1);
        p_upd=p(2:end,1);
end

%labels=M.stats(1).ea_stats.atlases.names(lidx);
edgecolor='w';
    jetlist=lines;

if exist('groupfiducial','var')
    edgecolor=jetlist(groupfiducial,:);
end

%jetlist(groups,:);
% plot areas:
f=figure('color','w','name',description,'Position',[100 100 550 550]);
g=gca;
if exist('color','var')
    scatter(g,X,Y,[],'o','MarkerEdgeColor',edgecolor,'MarkerFaceColor',color);
else
    hold on
    for group2=unique(groupfiducial)'
        edgecolor=jetlist(group2,:);
        scatter(g,X(groupfiducial==group2),Y(groupfiducial==group2),[],'o','MarkerEdgeColor',edgecolor,'MarkerFaceColor',jetlist(groupcolor(groupfiducial==group2),:));
    end
end

h=lsline;
set(h,'color','k');
%axis square

ax=gca;

title([description,' (R=',sprintf('%.3f',R_upd),', p=',sprintf('%.3f',p_upd),').'],'FontSize',16,'FontName','Helvetica');
xlabel(ea_underscore2space(labels{1}),'FontSize',16,'FontName','Helvetica');
ylabel(labels{2},'FontSize',16,'FontName','Helvetica');
axis square
% spacing=mean([nanvar(X(:,1)),nanvar(X(:,area+1))]);
% xlim([nanmin(X(:,1))-spacing,nanmax(X(:,2))+spacing]);
% ylim([nanmin(X(:,area+1))-spacing,nanmax(X(:,area+1))+spacing]);
