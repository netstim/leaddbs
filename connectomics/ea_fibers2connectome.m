function f=ea_fibers2connectome(fibs)
% small helper that will convert an nx4 matrix to lead-dbs connectome
% format 1.0

f.ea_fibformat='1.0';
f.fourindex=1;
f.fibers=fibs;

fidx=unique(fibs(:,4));
[h,ix]=ismember(fibs(:,4),fidx);
f.fibers(:,4)=ix; % make sure fourth column is 1:N format
cnt=1;
for fib=1:length(fidx)
    f.idx(cnt)=sum(ix==cnt);
    cnt=cnt+1;
end

