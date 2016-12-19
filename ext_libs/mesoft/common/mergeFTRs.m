function mergeFTRs(cellfile,cellrois,outname)


for k = 1:length(cellfile)
    ftr = ftrstruct_read(cellfile{k});
    id = 1:length(ftr.curveSegCell);
    if not(isempty(cellrois{k})),
        id = ftr.fiber{cellrois{k}}.curveID;        
    end;
    if k==1,
        mftr = ftr;
        mftr.curveSegCell = mftr.curveSegCell(id);
    else
        mftr.curveSegCell = {mftr.curveSegCell{:} ftr.curveSegCell{id}};
    end;
    
end;
mftr.fiber = {};
mftr.user = [];
mftr.connectCell = arrayfun(@(x) x, 1:length(mftr.curveSegCell),'uniformoutput',false)';
mftr.curveSegCell = mftr.curveSegCell';
ftrstruct_write(mftr,outname);
 
return;