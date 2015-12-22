function [fibdensRGB fibdens epdens] = ftr2FDmaps(ftr,sz,sel,M)


if not(isempty(sel)),
    idx = ftr.fiber{sel}.curveID;
else    
    idx = 1:length(ftr.curveSegCell);
end;
idx = idx(cellfun(@length,ftr.curveSegCell(idx))>1);
vox = ftr.vox;

nc = cellfun(@(x) ip(x,M*3),ftr.curveSegCell(idx),'UniformOutput',false);
nd = cellfun(@(x) id(x,M*3),ftr.curveSegCell(idx),'UniformOutput',false);

for k = 1:length(nc),
    dir = nd{k};
    dir(:,1) = dir(:,1)*vox(1);
    dir(:,2) = dir(:,2)*vox(2);
    dir(:,3) = dir(:,3)*vox(3);
    arclen = sqrt(sum(dir.^2,2));
    ncs{k} = dir;
    ncal{k} = arclen;
end;

pts  = cat(1,nc{:})-1;
dir = cat(1,ncs{:});
arclen = cat(1,ncal{:});
perm = [2 1 3];
for k = 1:3,
    fibdensRGB(:,:,:,perm(k)) = AccumulateBilinWeighted(single(sz*M),single(pts'*M),single(abs(dir(:,k))));
end;

fibdens = AccumulateBilinWeighted(single(sz*M),single(pts'*M),single(arclen)); 


nc = cellfun(@ep,ftr.curveSegCell(idx),'UniformOutput',false);
pts  = cat(1,nc{:})-1;
epdens = AccumulateBilin(single(sz*M),single(pts'*M)); 



return;


function y = ip(x,N)
dx = x(2:end,:) - x(1:end-1,:);
x = x(1:end-1,:);
y = ones(N,size(x,1),size(x,2));
for k = 1:N,
    y(k,:,:) = x + k/N*dx ;
end;
y = reshape(y,[N*size(x,1) size(x,2)]);

function y = id(x,N)
dx = x(2:end,:) - x(1:end-1,:);
x = x(1:end-1,:);
y = ones(N,size(x,1),size(x,2));
for k = 1:N,
    y(k,:,:) = dx ;
end;
y = reshape(y,[N*size(x,1) size(x,2)]);



function y = ep(x)
y = [x(1,:) ; x(end,:)];
