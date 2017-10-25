function tmct=ea_tonemap_ct(ct)
tmct=ct;
% if min(tmct(:))<-5000
%     disp('CT seems to have incorrect Hounsfield units, attempting datadriven normalization.')
%     tmct(:)=ea_nanzscore(tmct(:));
%     try
%     normvals=ea_normal(tmct(1:1000:end)');
%     f=fit(tmct(1:1000:end)',normvals,'poly2');
%     tmct(:)=f(tmct(:));
%     end
%     %figure, hist(normvals,1000)
%     return
% end
% brain window: center = 40, width = 80
tmct(ct>0 & ct<80)=tmct(ct>0 & ct<80)/80;

% bone window: center = 300, width = 2000
tmct(ct>=80 & ct<1300)=(tmct(ct>=80 & ct<1300)-80)/1220;

% saturate above and below levels:
tmct(ct>=1300)=1;
tmct(ct<0)=0;
