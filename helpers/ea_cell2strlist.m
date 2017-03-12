function strlist=ea_cell2strlist(incell)
strlist=[];
for c=1:length(incell)
    strlist=[strlist,incell{c},', '];
end
strlist=strlist(1:end-2);
