function [h,p]=ea_plothistperm(title,similarities,idlabels,ids,cols,onesided)
if ~exist('cols','var') || isempty(cols)
   cols=ones(length(ids),1);
end
if ~exist('onesided','var')
    onesided=0;
end

% remove true sims from similarities:
todel=[];
R0similarities=similarities;
for id=1:length(ids)
   todel=[todel;ids{id}];
end
R0similarities(todel)=[];


colorc=ones(length(R0similarities),1);

gr(1)=gramm('x',R0similarities,'color',colorc);
gr(1).set_color_options('chroma',0,'lightness',50);
gr(1).stat_bin('geom','line','fill','all');

gr(1).set_title(['> ',title]);
h=figure('Position',[100 100 800 600]);

gr.draw();
hold on

set(0,'CurrentFigure',h);
set(h,'CurrentAxes',gr.facet_axes_handles);

ccode=lines;

% plot lines
for g=1:length(idlabels)
    thissims=similarities(ids{g});

    plot([mean(thissims),mean(thissims)],[0,10000],'Color',ccode(cols(g),:),'LineWidth',2);

end
ax=gr(1).facet_axes_handles;

% plot textboxes
for g=1:length(idlabels)
    thissims=similarities(ids{g});
    if onesided
        ssims=sort((R0similarities),'descend');
        ssims=ssims>(mean(thissims));
    else
        ssims=sort(abs(R0similarities),'descend');
        ssims=ssims>abs(mean(thissims)); 
    end
    p=sum(ssims)/numel(ssims);
    if mean(thissims)>0
        text(double(mean(thissims)),(ax.YLim(2))-(g*(0.07*ax.YLim(2))),[' \leftarrow ',ea_underscore2space(idlabels{g}),' [p = ',sprintf('%.4f',p),']'],'Color',ccode(cols(g),:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','left','BackgroundColor','w');
    else
        text(double(mean(thissims)),(ax.YLim(2))-(g*(0.07*ax.YLim(2))),['',ea_underscore2space(idlabels{g}),' [p = ',sprintf('%.4f',p),'] \rightarrow '],'Color',ccode(cols(g),:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','right','BackgroundColor','w');
    end
end
