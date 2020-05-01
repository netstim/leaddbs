function ea_ttest(x,y,description,labels)
% this simple function is a small wrapper for a t-test figure and has been
% written to go in line with the ea_corrplot wrapper.

disp(description);

[h,p,ci,stats]=ttest2(x,y)

ygramm=[x;y]; % concat for gramm
xgramm=labels([ones(length(x),1);repmat(2,length(y),1)]);
g(1)=gramm('x',ones(length(ygramm),1),'y',ygramm,'color',xgramm);
g(1).stat_boxplot();
g(1).set_title(description);

%boxplot([stats.tstat+ci(1),stats.tstat,stats.tstat+ci(2)])
g.set_names('x','Group','y','Values','color','Group');

figure('Position',[27   449   312   571]);
g.draw();
