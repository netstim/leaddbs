function h=ea_checkgenesim(map,smooth,geneidx,cols)
prefs=ea_prefs;
if ~iscell(geneidx)
    geneidx={geneidx};
end

load([prefs.genetics.dbdir,'geneinfo.mat'],'geneinfo');
for g=1:length(geneidx)
    id{g}=[];
    tableheaders=fieldnames(geneinfo);
    for header=[2,4,5]
        id{g}=ismember(geneinfo.(tableheaders{header}),geneidx{g});
        if any(id{g})
            break
        end
    end
    if ~any(id{g})
        ea_error(['Could not find gene: ',geneidx{g},'.']);
    end
    id{g}=find(id{g});
end

if iscell(map)
    for m=1:length(map)
        simsm(m,:)=ea_probeexpressions(map{m},[prefs.genetics.dbdir,'genedb.mat'],[],smooth(m));
    end
    sims=sum(simsm,1);
%     for m=1:length(map)
%         if ~exist('sims','var')
%             sims=simsm(1,:);
%         else
%             sims=sims.*(simsm(m,:));
%         end
%     end
    [mapfn]='Multimap';
else
    sims=ea_probeexpressions(map,[prefs.genetics.dbdir,'genedb.mat'],[],smooth);
    [~,mapfn]=fileparts(map);
end

if ~exist('cols','var')
    ea_plothistperm(mapfn,sims,geneidx,id);
else
    ea_plothistperm(mapfn,sims,geneidx,id,cols);
end
% colorc=ones(length(sims),1);
% 
% gr(1)=gramm('x',sims,'color',colorc);
% gr(1).set_color_options('chroma',0,'lightness',50);
% gr(1).stat_bin('geom','line','fill','all');
% 
% gr(1).set_title(['> ',mapfn]);
% h=figure('Position',[100 100 800 600]);
% 
% gr.draw();
% hold on
% 
% set(0,'CurrentFigure',h);
% set(h,'CurrentAxes',gr.facet_axes_handles);
% 
% ccode=lines;
% 
% % plot lines
% for g=1:length(geneidx)
%     thissims=sims(id{g});
% 
%     plot([mean(thissims),mean(thissims)],[0,10000],'Color',ccode(cols(g),:),'LineWidth',2);
%     
% end
% 
% ax=gr(1).facet_axes_handles;
% 
% 
% % plot textboxes
% for g=1:length(geneidx)
%     thissims=sims(id{g});
% 
%     ssims=sort(abs(sims),'descend');
%     
%     ssims=ssims>abs(mean(thissims));
%     p=sum(ssims)/numel(ssims);
%     if mean(thissims)>0
%         text(double(mean(thissims)),(ax.YLim(2))-(g*(0.07*ax.YLim(2))),[' \leftarrow ',geneidx{g},' [p = ',sprintf('%.4f',p),']'],'Color',ccode(cols(g),:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','left','BackgroundColor','w');
%     else
%         text(double(mean(thissims)),(ax.YLim(2))-(g*(0.07*ax.YLim(2))),['',geneidx{g},' [p = ',sprintf('%.4f',p),'] \rightarrow '],'Color',ccode(cols(g),:),'FontSize',14,'FontWeight','bold','HorizontalAlignment','right','BackgroundColor','w');
%     end
% end
