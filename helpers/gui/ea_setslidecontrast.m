function ea_setslidecontrast(~,~,contrastoffset,posneg,resultfig,handles)
sc=getappdata(resultfig,'slidecontrast');
if isempty(sc)
    sc.c=0;
    sc.o=0;
end
switch contrastoffset
    case 'c' % contrast
        sc.c=sc.c+posneg;
    case 'o' % offset
        sc.o=sc.o+posneg*0.2;
end
setappdata(resultfig,'slidecontrast',sc);


togglestates = getappdata(resultfig,'togglestates');
togglestates.refreshview = 1;
setappdata(resultfig,'togglestates',togglestates);

if ~exist('handles','var')
   handles=resultfig; 
end

ea_refreshresultfig(handles,1);