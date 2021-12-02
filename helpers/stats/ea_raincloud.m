function [h,T,p,g]=ea_raincloud(X,Y,labels,corrtype,group1,group2,pperm,colors,markers)

% X is the dividing (binary) var while Y is the real deal
[hyp,p,ci,stats]=ttest2(Y(X==0),Y(X==1));
T=stats.tstat;

h=figure('Name',labels{1},'Color','w','NumberTitle','off');

sp1=subplot(2,1,1);
subtitle(labels{2});
col=ea_color_wes('lifeaquatic');
g=ea_raincloud_plot(Y(X==0),'color',col(1,:),'box_on',1);
a1=gca;
set(a1,'ytick',[])
a1.YLabel.String='Control';
a1.YLabel.FontSize=14;

a1.Box='off';

text(0.8,0.8,['T = ',sprintf('%0.2f',T),'; p = ',sprintf('%0.2f',p)],'FontWeight','bold','FontSize',14,'HorizontalAlignment','right','Units','normalized');

sp1=subplot(2,1,2);
subtitle(labels{2});
g=ea_raincloud_plot(Y(X==1),'color',col(3,:),'box_on',1);
a2=gca;
set(a2,'ytick',[])
a2.YLabel.String=labels{2};
a2.YLabel.FontSize=14;
a2.Box='off';

