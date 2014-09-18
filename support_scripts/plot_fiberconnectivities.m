function plot_fiberconnectivities(M)


clc
close all

% gather fiber connectivities
for side=1:2
    switch side
        case 1
            lidx=[2:2:90];
        case 2
            lidx=[1:2:90];
    end
    
    
    for pt=1:length(M.patient.list)
        cb=[mean(M.stats(pt).ea_stats.stimulation(end).ft(side).nfibercounts{1}(91:end))]';
        fc_c{side}(pt,:)=[M.stats(pt).ea_stats.stimulation(end).ft(side).nfibercounts{1}(lidx);cb]';
    end
    fcs=fc_c{side};


% process fcs.
fcmean=max(fcs);
scalef=max(fcmean);

[fcmean,idx]=sort(fcmean,2,'descend');


fcs=fcs(:,idx);

h=figure;
set(h,'Position',[200,200,1000,600]);

set(h,'Color','w');
for pt=1:length(M.patient.list)
    cvarnnan=M.clinical.vars{1};
    cdat=repmat((M.clinical.vars{1}(pt)-min(M.clinical.vars{1}))/(max(M.clinical.vars{1})-min(M.clinical.vars{1})),1,3);
    p(pt)=plot(fcs(pt,:)','Color',cdat,'LineSmoothing','on','LineWidth',2);

hold('on')
end



[R,p]=corrcoef([M.clinical.vars{1},fcs],'rows','pairwise');
p_upd=p(2:end,1);
labels=[M.stats(1).ea_stats.stimulation.ft(1).labels{1}(lidx);{'Cerebellum_B'}];
labels=labels(idx);



for fc=1:length(fcmean);
    if p_upd(fc)<0.05
        cwsig='r';
        append='';
    elseif p_upd(fc)<0.001
        cwsig='r';
        append='*';
    else
        cwsig='k';
        append='';
    end
    t(fc)=text(fc,fcmean(fc)+scalef/20,[sub2space(labels{fc}(1:end-2)),append],'FontName','Helvetica','FontSize',12,'color',cwsig);
    plot([fc,fc],[fcmean(fc),fcmean(fc)+scalef/20],'Color',[0.7,0.7,0.7]);
    set(t(fc),'rotation', 55);
    if abs(diff([fcs(1,fc),fcs(2,fc)]))>0.5*max([fcs(1,fc),fcs(2,fc)]) && mean([fcs(1,fc),fcs(2,fc)])>10
        set(t(fc),'FontWeight','bold');
        set(t(fc),'Color','r');
        plot([fc,fc],[fcs(1,fc),fcs(2,fc)],'r*-');
    end
end
%legend(p,'K0','K1','K2','K3','K8','K9','K10','K11','Location','best');

axis([0 fc+10 0 max(fcmean)+scalef/3]);
ti = get(gca,'TightInset');
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)])
box off
set(gca,'XTick',[])
saveas(h,['fiberconnectivities',num2str(side),'.png']);


%take_screenshot;
end


function str=sub2space(str) % replaces subscores with spaces
str(str=='_')=' ';