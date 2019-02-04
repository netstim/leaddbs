function [h,pv,Rsquared,F,mdl]=ea_glmplot(X,y,labels,distribution)

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


mdl=fitglm(X,y,'distribution',distribution);

yhat=predict(mdl,X);
g=gramm('x',yhat,'y',y); % data needs to be put in "reversed" for gramm.
g.geom_point();
g.stat_glm('distribution',distribution,'fullrange','false','fullrange','false');



%g.geom_abline();
%g.stat_cornerhist('edges',[ea_nanmean(X)-1*ea_nanstd(X):0.2:ea_nanmean(X)+1*ea_nanstd(X)],'aspect',0.6,'location',max(X));
pstr='p';
[pv,F]=coefTest(mdl);
Rsquared=mdl.Rsquared.Ordinary;
g.set_title([labels{1},' [R2 = ',sprintf('%.2f',Rsquared),'; ',pstr,' = ',sprintf('%.3f',pv),'; F-stat = ',sprintf('%.2f',F),']'],'FontSize',20);
g.set_names('x',labels{2},'y',labels{3});
g.set_text_options('base_size',22);
g.no_legend();
h=figure('Position',[100 100 550 550]);
g.draw();
