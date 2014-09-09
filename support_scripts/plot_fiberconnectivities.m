function plot_fiberconnectivities(M)

keyboard
clc
close all
load('ea_stats.mat');

% gather fiber connectivities
for pt=1:length(M.patient.list)

fcs(pt,:)=[ea_stats.ft(1).fibercounts{1}(2:2:90)'
    ea_stats.ft(2).fibercounts{1}(2:2:90)'
    ea_stats.ft(3).fibercounts{1}(2:2:90)'
    ea_stats.ft(4).fibercounts{1}(2:2:90)'

    ea_stats.ft(5).fibercounts{1}(1:2:90)'
    ea_stats.ft(6).fibercounts{1}(1:2:90)'
    ea_stats.ft(7).fibercounts{1}(1:2:90)'
    ea_stats.ft(8).fibercounts{1}(1:2:90)'
];


cb=[mean(ea_stats.ft(1).fibercounts{1}(91:end))
    mean(ea_stats.ft(2).fibercounts{1}(91:end))
    mean(ea_stats.ft(3).fibercounts{1}(91:end))
    mean(ea_stats.ft(4).fibercounts{1}(91:end))
    mean(ea_stats.ft(5).fibercounts{1}(91:end))
    mean(ea_stats.ft(6).fibercounts{1}(91:end))
    mean(ea_stats.ft(7).fibercounts{1}(91:end))
    mean(ea_stats.ft(8).fibercounts{1}(91:end))];

fcs=[fcs,cb]; % add cerebellum as one entry.

end

% process fcs.
fcmean=max(fcs);
scalef=max(fcmean);

[fcmean,idx]=sort(fcmean,2,'descend');


fcs=fcs(:,idx);

h=figure;
set(h,'Position',[200,200,1000,600]);

set(h,'Color','w');
p(1)=plot(fcs(1,:)','Color',[0.1,0.1,1],'LineSmoothing','on','LineWidth',2);
hold('on')
p(2)=plot(fcs(2,:)','Color',[0.1,0.1,0.9],'LineSmoothing','on','LineWidth',2);
p(3)=plot(fcs(3,:)','Color',[0.1,0.1,0.8],'LineSmoothing','on','LineWidth',2);
p(4)=plot(fcs(4,:)','Color',[0.1,0.1,0.7],'LineSmoothing','on','LineWidth',2);
p(5)=plot(fcs(5,:)','Color',[1,0.1,0.1],'LineSmoothing','on','LineWidth',2);
p(6)=plot(fcs(6,:)','Color',[0.9,0.1,0.1],'LineSmoothing','on','LineWidth',2);
p(7)=plot(fcs(7,:)','Color',[0.8,0.1,0.1],'LineSmoothing','on','LineWidth',2);
p(8)=plot(fcs(8,:)','Color',[0.7,0.1,0.1],'LineSmoothing','on','LineWidth',2);


labels=[ea_stats.ft(1).labels{1}(1:2:90);{'Cerebellum_B'}];
labels=labels(idx);



for fc=1:length(fcmean);
    t(fc)=text(fc,fcmean(fc)+scalef/20,sub2space(labels{fc}(1:end-2)),'FontName','Helvetica','FontSize',12);
    plot([fc,fc],[fcmean(fc),fcmean(fc)+scalef/20],'Color',[0.7,0.7,0.7]);
    set(t(fc),'rotation', 55);
    if abs(diff([fcs(1,fc),fcs(2,fc)]))>0.5*max([fcs(1,fc),fcs(2,fc)]) && mean([fcs(1,fc),fcs(2,fc)])>10
        set(t(fc),'FontWeight','bold');
        set(t(fc),'Color','r');
        plot([fc,fc],[fcs(1,fc),fcs(2,fc)],'r*-');
    end
end
legend(p,'K0','K1','K2','K3','K8','K9','K10','K11','Location','best');

axis([0 fc+10 0 max(fcmean)+scalef/3]);
ti = get(gca,'TightInset');
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)])
box off
set(gca,'XTick',[])
saveas(h,'fiberconnectivities.png');


take_screenshot;


function str=sub2space(str) % replaces subscores with spaces
str(str=='_')=' ';