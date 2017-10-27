function tmct=ea_tonemap_ct(ct)
prefs=ea_prefs;
tmct=ct;
switch prefs.tonemap
    case 'heuristic'

tmct(ct>0 & ct<80)=tmct(ct>0 & ct<80)/80;

% bone window: center = 300, width = 2000
tmct(ct>=80 & ct<1300)=(tmct(ct>=80 & ct<1300)-80)/1220;

% saturate above and below levels:
tmct(ct>=1300)=1;
tmct(ct<0)=0;
    case 'albada'
        tmct(:)=ea_normal(ct(:),1,0,' ',0,1,'TRUE');     
end
