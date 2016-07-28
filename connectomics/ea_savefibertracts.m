function ea_savefibertracts(cfile,fibers,idx,voxmm,mat)


[pth,fn,~]=fileparts(cfile);
ftr.fourindex=1;
ftr.ea_fibformat='1.1';
ftr.fibers=fibers;
ftr.idx=idx;
ftr.voxmm=voxmm;
if nargin == 5
    ftr.mat=mat;
end
display(sprintf('\nSaving fibers: %s.mat...',fn));
save(fullfile(pth,[fn,'.mat']),'-struct','ftr','-v7.3');
disp('Done.');
