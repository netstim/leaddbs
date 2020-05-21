function ea_savefibertracts(outputFTR,fibers,idx,voxmm,mat,vals)

ftr.fourindex=1;
ftr.ea_fibformat='1.1';
ftr.fibers=fibers;
ftr.idx=idx;
ftr.voxmm=voxmm;

if exist('mat','var')
    if ischar(mat) % Path to the reference image provided
        ftr.mat = ea_get_affine(mat);
    else
        ftr.mat = mat;
    end
end

if exist('vals','var')
    ftr.vals=vals;
else
    ftr.vals=ones(size(ftr.idx));
end

fprintf('\nSaving fibers to %s...\n', outputFTR);
save(outputFTR,'-struct','ftr','-v7.3');
disp('Done.');
