function ea_ttest(x,y,description,labels)
% this simple function is a small wrapper for a t-test figure and has been
% written to go in line with the ea_corrplot wrapper.

disp(description);

[h,p,ci,stats]=ttest2(x,y)

figure

set(gcf,'name',description);
set(gcf,'NumberTitle','off');

boxplot([stats.tstat+ci(1),stats.tstat,stats.tstat+ci(2)])

title(['t-value: ',num2str(stats.tstat),', ','p-value: ',num2str(p),'.']);

