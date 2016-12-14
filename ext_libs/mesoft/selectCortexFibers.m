function ftr = selectCortexFibers(finame,gm,oversamp,Nsz)

gm = nifti(gm);
gm = double(gm.dat);

ftr = ftrstruct_read(finame);

fibs = ftr.curveSegCell;
terms = cellfun(@(x) [x(1,:); x(end,:)]-1,fibs,'UniformOutput',false);
terms = cat(1,terms{:})'*oversamp;

idxbol  = SelectCorticalFibers(single(gm),double(Nsz),single(terms));
ftr.connectCell = ftr.connectCell(idxbol>0);
ftr.curveSegCell = ftr.curveSegCell(idxbol>0);
ftr.curveD = ftr.curveD(idxbol>0);

