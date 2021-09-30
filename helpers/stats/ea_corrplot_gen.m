function [R_upd,p_upd,R,p,f]=ea_corrplot_gen(varargin)
% this simple function is a small wrapper for a corrplot figure.
% [R_upd,p_upd,R,p,f]=ea_corrplot(X,description,labels,handles,color/groups,corrtype)

X=varargin{1};
description=varargin{2};
labels=varargin{3};
if nargin>3
    handles=varargin{4};
end
if nargin>4
    if ~isempty(varargin{5})
        if ischar(varargin{5}) || isequal([1 3],size(varargin{5}))
            color=varargin{5};
        else
            groups=varargin{5};
        end
    else
        groups=ones(length(X),1);
    end
else
    groups=ones(length(X),1);
end

corrtype='Pearson';
try corrtype=varargin{6}; end

dim=size(X,2);

disp(description);


nnans=isnan(sum(X,2));
X(nnans,:)=[];
try groups(nnans)=[]; end

switch corrtype
    case {'permutation_spearman','permutation'}
        [R,p]=ea_permcorr(X(:,1),X(:,2),'spearman');
        R_upd=R; p_upd=p;
    case 'permutation_pearson'
        [R,p]=ea_permcorr(X(:,1),X(:,2),'pearson');
        R_upd=R; p_upd=p;
    otherwise
[R,p]=corr(X,'rows','pairwise','type',corrtype);
R_upd=R(2:end,1);
p_upd=p(2:end,1);
end

% support for nans
%labels=M.stats(1).ea_stats.atlases.names(lidx);
edgecolor='w';
if size(groups,2)>1
    edgecolor=jetlist(groups(:,2),:);
    groups=groups(:,1);
end

%jetlist(groups,:);
for area=1:length(R_upd)
    %% plot areas:
    f=figure('color','w','name',description,'Position',[100 100 550 550]);
    jetlist=lines;
    g=gca;
    if exist('color','var')
        scatter(g,X(:,1),X(:,area+1),[],'o','MarkerEdgeColor',edgecolor,'MarkerFaceColor',color);
    else
        try
            scatter(g,X(:,1),X(:,area+1),[],'o','MarkerEdgeColor',edgecolors,'MarkerFaceColor',jetlist(groups,:));
        catch
            scatter(g,X(:,1),X(:,area+1),[],jetlist(groups,:),'filled');
        end
    end

    h=lsline;
    set(h,'color','k');
    %axis square

    ax=gca;

    [~,fn]=fileparts(labels{area+1});
    if length(fn)>4
        if strcmp(fn(end-3:end),'.nii')
            [~,fn]=fileparts(fn);
        end
    end

    title([description,' (R=',sprintf('%.3f',R_upd(area)),', p=',sprintf('%.3f',p_upd(area)),').'],'FontSize',16,'FontName','Helvetica');
    xlabel(ea_underscore2space(labels{1}),'FontSize',16,'FontName','Helvetica');
    ylabel(labels{2},'FontSize',16,'FontName','Helvetica');
    axis square
    % spacing=mean([nanvar(X(:,1)),nanvar(X(:,area+1))]);
    % xlim([nanmin(X(:,1))-spacing,nanmax(X(:,2))+spacing]);
    % ylim([nanmin(X(:,area+1))-spacing,nanmax(X(:,area+1))+spacing]);
    if nargin>3
        if ~isempty(varargin{4})
            odir=handles.groupdir_choosebox.String;
            ofname=[odir,description,'_',fn,'_',labels{1},'.png'];
            ea_screenshot(ofname);
        end
    end
end
