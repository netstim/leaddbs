function ea_savefibertracts(cfile,fibers,idx,voxmm,mat,vals)

[pth,fn,~]=fileparts(cfile);
ftr.fourindex=1;
ftr.ea_fibformat='1.1';
ftr.fibers=fibers;
ftr.idx=idx;
ftr.voxmm=voxmm;
if exist('mat','var')
    ftr.mat=mat;
end
if exist('vals','var')
    ftr.vals=vals;
else
    ftr.vals=ones(size(ftr.idx));
end
fprintf('\nSaving fibers: %s.mat...\n',fn);
save(fullfile(pth,[fn,'.mat']),'-struct','ftr','-v7.3');
disp('Done.');
