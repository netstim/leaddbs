function [res ftr] = mesoftr2FDmaps(ftr,M)

if isstr(ftr),
    ftr = load(ftr);
end;

if not(exist('M')),
    M = 1;
end;

sz = size(ftr.user.S2(:,:,:,1));
idx = 1:length(ftr.curveSegCell);
idx = idx(cellfun(@length,ftr.curveSegCell(idx))>1);

if isempty(idx)
    res = [];
    return;
end;

vox = ftr.vox;

t = 3;

nc = cellfun(@(x) ip(x,M*t),ftr.curveSegCell(idx),'UniformOutput',false);
nd = cellfun(@(x) id(x,M*t),ftr.curveSegCell(idx),'UniformOutput',false);
if isfield(ftr,'curveD'),
    pa = cellfun(@(x) ip(x',M*t),ftr.curveD(idx),'UniformOutput',false);
end;


for k = 1:length(nc),
    dir = nd{k};
    dir(:,1) = dir(:,1)*vox(1);
    dir(:,2) = dir(:,2)*vox(2);
    dir(:,3) = dir(:,3)*vox(3);
    arclen = sqrt(sum(dir.^2,2));
    ncs{k} = dir;
    ncal{k} = arclen;
end;

pts  = cat(1,nc{:})-1.5;
if isfield(ftr,'curveD'),
    pa = cat(1,pa{:});
end;
dir = cat(1,ncs{:});
arclen = cat(1,ncal{:});

perm = [2 1 3];
for k = 1:3,
    fibdensRGB(:,:,:,perm(k)) = AccumulateTrilinWeighted(single(sz*M),single(pts'*M),single(abs(dir(:,k))));
end;
fibdensRGB = fibdensRGB /max(fibdensRGB(:));

fibdens = AccumulateTrilinWeighted(single(sz*M),single(pts'*M),single(arclen));
if isfield(ftr,'curveD'),
    for k = 1:size(pa,2),
        PaMap{k} = AccumulateTrilinWeighted(single(sz*M),single(pts'*M),single(pa(:,k).*arclen)')./(eps+fibdens);
    end;
end;



nc = cellfun(@ep,ftr.curveSegCell(idx),'UniformOutput',false);
pts  = cat(1,nc{:})-1;
epdens = AccumulateTrilin(single(sz*M),single(pts'*M));

res.fdrgb = fibdensRGB;
res.fd = fibdens;
res.ep = epdens;
if isfield(ftr,'curveD'),
    maps = computeParameterMaps(ftr,M);
    res.Daxon = PaMap{1};
    res.Dextrapara = PaMap{2};
    res.Dextraorth = PaMap{3};
    res.turt = PaMap{2}./(PaMap{3}+eps);
    res.vfi = maps.vfi;
    res.vfsw = maps.vfsw;
    res.vfe = maps.vfe;
    res.S0 = maps.S0;
    res.segcount = maps.segcount;

    gauss = fspecial('gaussian',[5 1],1); [gx gy gz] = ndgrid(gauss); gauss = gx.*gy.*gz;
    res.axondiameter = maps.vfi ./ (imfilter(maps.segcount,gauss)./(eps+ imfilter(double(maps.segcount>0),gauss)));



end;


Q = diag([1/M 1/M 1/M 1]);
mr.edges = ftr.hMatrix*Q;
mr.vox = ftr.vox/M;

res.mrProp = mr;

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






