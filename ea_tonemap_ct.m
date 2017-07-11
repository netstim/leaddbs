function tmct=ea_tonemap_ct(ct)
tmct=ct;

% brain window: center = 40, width = 80
tmct(ct>0 & ct<80)=tmct(ct>0 & ct<80)/80;

% bone window: center = 300, width = 2000
tmct(ct>80 & ct<1300)=(tmct(ct>80 & ct<1300)-80)/1220;

% saturate above and below levels:
tmct(ct>=1300)=1;
tmct(ct<0)=0;