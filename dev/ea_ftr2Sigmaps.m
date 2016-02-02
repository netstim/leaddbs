function [recompSig vfi] = ftr2FDmaps(ftr,bTensor)
% function by M. Reisert to estimate diffusion signal based on fiber tracts
% that have been estimated using the mesoFT tracking system.

sel = [];
M = 1;
bTensor = bTensor/1000;

if not(isempty(sel)),
    idx = ftr.fiber{sel}.curveID;
else    
    idx = 1:length(ftr.normalized_fibers_vox);
end;
idx = idx(cellfun(@length,ftr.normalized_fibers_vox(idx))>1);


   osamp =  (ftr.trackParam.params.p_wid);
t = 5;

nc = cellfun(@(x) ip(x,M*t),ftr.normalized_fibers_vox(idx),'UniformOutput',false);
nd = cellfun(@(x) id(x,M*t),ftr.normalized_fibers_vox(idx),'UniformOutput',false);
pa = cellfun(@(x) ip(x',M*t),ftr.curveD(idx),'UniformOutput',false);
pts  = ((cat(1,nc{:})-0.5)*osamp);
dirs  = (cat(1,nd{:}));
pa = cat(1,pa{:});

 pts = ftr.user.P(1:3,:);
 pts = diag(osamp./ftr.vox)*pts;
 pts = floor(pts');
 dirs = ftr.user.P(4:6,:);
 dirs = dirs';
 pa = ftr.user.P(11:15,:);
 pa = pa';



vfi = interp3(ftr.user.vf(:,:,:,1),floor(pts(:,2))+1,floor(pts(:,1))+1,floor(pts(:,3))+1);
szvf = size( ftr.user.vf);
recompSig = zeros([szvf(1:3)*M,size(bTensor,3)]);
for k = 1:size(bTensor,3),
    q = dirs(:,1).^2*bTensor(1,1,k) + dirs(:,2).^2*bTensor(2,2,k) +  dirs(:,3).^2*bTensor(3,3,k) + 2*(dirs(:,1).*dirs(:,2).*bTensor(1,2,k) + dirs(:,3).*dirs(:,2).*bTensor(3,2,k) + dirs(:,1).*dirs(:,3).*bTensor(1,3,k));
    tr = bTensor(1,1,k) + bTensor(2,2,k) + bTensor(3,3,k);
    S = (pa(:,4).*(vfi.*exp(-pa(:,1).*q) + (1-vfi).*exp(-(pa(:,2).*q + (tr-q).*pa(:,3)))));
    recompSig(:,:,:,k) = AccumulateBilinWeighted(single(szvf*M),single(pts'*M),single(S));            
end;

tr = squeeze(bTensor(1,1,:) + bTensor(2,2,:) + bTensor(3,3,:));
b0 = tr==0;
b0avg = mean(recompSig(:,:,:,b0),4);
recompSig = imfilter(recompSig,ones(osamp,osamp,osamp));
recompSig = recompSig(1:osamp:end,1:osamp:end,1:osamp:end,:);
% for k = 1:size(bTensor,3),
%     recompSig(:,:,:,k) = recompSig(:,:,:,k) ./ b0avg;
% end;




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
