function [fibers,idx]=ea_loadfibertracts(cfile)

fibinfo=load(cfile);
if ~isfield(fibinfo,'ea_fibformat')
fibinfo=ea_convertfibs2newformat(fibinfo,cfile);
end

fibers=fibinfo.fibers;
idx=fibinfo.idx;





function ftr=ea_convertfibs2newformat(fibinfo,cfile)
    fn=fieldnames(fibinfo);

if isfield(fn,'normalized_fibers_mm')
    fibers=fibinfo.normalized_fibers_mm;
else
    fibers=eval(['fibinfo.',fn{1},';']);
end

c=size(fibers);
if c(1)<c(2)
    fibers=fibers';
end

ea_dispercent(0,'Converting fibers');
[idx,~]=cellfun(@size,fibers);

fibers=cell2mat(fibers);
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'
    ea_dispercent(cnt/length(idx));
    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
ea_dispercent(1,'end');

fibers=[fibers,idxv];

[pth,fn,ext]=fileparts(cfile);
ftr.fourindex=1;
ftr.ea_fibformat='1.0';
ftr.fibers=fibers;
ftr.idx=idx;
disp('Saving fibers in new format');
save(fullfile(pth,[fn,'.mat']),'-struct','ftr','-v7.3');
disp('Done.');