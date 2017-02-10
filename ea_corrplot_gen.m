function [R_upd,p_upd,R,p,f]=ea_corrplot(varargin)
% this simple function is a small wrapper for a corrplot figure.


X=varargin{1};
description=varargin{2};
labels=varargin{3};
if nargin>3
    handles=varargin{4};
end
if nargin>4
    if ischar(varargin{5}) || isequal([1 3],size(varargin{5}))
        color=varargin{5};
    else
        groups=varargin{5};
    end
else
    groups=ones(length(X),1);
end

dim=size(X,2);

disp(description);




[R,p]=corr(X,'rows','pairwise','type','Pearson');
R_upd=R(2:end,1);
p_upd=p(2:end,1);

%labels=M.stats(1).ea_stats.atlases.names(lidx);


jetlist=lines;
%jetlist(groups,:);
for area=1:length(R_upd)
    %% plot areas:
    f=figure('color','w','name',description);
    if exist('color','var')
        scatter(X(:,1),X(:,area+1),[],color,'filled');
    else
        
        scatter(X(:,1),X(:,area+1),[],jetlist(groups,:),'filled');

        
    end
    
    h=lsline;
    set(h,'color','k');
    axis square
    
    ax=gca;
% ax.YTick=[-20:20:90];
% ax.YLim=[-15,90];
% ax.XTick=[-20:10:90];
% ax.XLim=[25,65];
[~,fn]=fileparts(labels{area+1});
    if strcmp(fn(end-3:end),'.nii')
        [~,fn]=fileparts(fn);
    end
    title([description,' (R=',sprintf('%.3f',R_upd(area)),', p=',sprintf('%.3f',p_upd(area)),').'],'FontSize',16,'FontName','Helvetica');
    xlabel(sub2space(labels{1}),'FontSize',16,'FontName','Helvetica');
    ylabel(labels{2},'FontSize',16,'FontName','Helvetica');
 %   spacing=mean([nanvar(X(:,1)),nanvar(X(:,area+1))]);
 %   xlim([nanmin(X(:,1))-spacing,nanmax(X(:,2))+spacing]);
 %   ylim([nanmin(X(:,area+1))-spacing,nanmax(X(:,area+1))+spacing]);
    if nargin==4
        if ~isempty(varargin{4})
        odir=get(handles.groupdir_choosebox,'String');
        ofname=[odir,description,'_',fn,'_',labels{1},'.png'];
        ea_screenshot(ofname);
        end
    end
end



function str=sub2space(str) % replaces subscores with spaces
str(str=='_')=' ';