function [signal themask] = genPhan(ten,nz)

fac = 2;
n = 12*fac;
d = 3;

[X Y] = ndgrid(1:n,1:n);
R = sqrt((X-n/2).^2 + (Y-n/2).^2);

cnt = 1;


Dcsf = 2;
Dp = 0.5;
vf = 0.7;
vsw = 0.0;

t = 0.1;
  
Di = 2.4;
Da = 1.7;
Dp = 0.6;


P(cnt).Di = Di;
P(cnt).Da = Da;
P(cnt).Dp = Dp;
P(cnt).vf = vf;
P(cnt).vsw = vsw;
mask{cnt} = R<n*0.4 & R >n*0.25;
dirs{cnt} = cat(3,(Y-n/2)./(eps+R),-(X-n/2)./(eps+R));
 cnt = cnt + 1;

P(cnt).Di = Di;
P(cnt).Da = Da; 
P(cnt).Dp = Dp;
P(cnt).vf = vf; %0.2+0.6*X/n;
P(cnt).vsw = vsw;
v = [1 -0.1]; v = v/norm(v);
%v = [-0.6 1]; v = v/norm(v);
p = ((X-n/2)*v(1) + (Y-n/2)*v(2));
mask{cnt} = p < n*t & p > -n*t;
dirs{cnt} = cat(3,-v(2)*ones(size(X)),v(1)*ones(size(X)));
 cnt = cnt + 1;

P(cnt).Di = Di;
P(cnt).Da = Da; 
P(cnt).Dp = Dp;
P(cnt).vf = vf; %0.2+0.6*X/n;
P(cnt).vsw = vsw;
%v = [1 -0.1]; v = v/norm(v);
v = [-0.6 1]; v = v/norm(v);
p = ((X-n/2)*v(1) + (Y-n/2)*v(2));
mask{cnt} = p < n*t & p > -n*t;
dirs{cnt} = cat(3,-v(2)*ones(size(X)),v(1)*ones(size(X)));
cnt = cnt + 1;


% P(cnt).Di = P(cnt-1).Di;
% P(cnt).Da = P(cnt-1).Da;
% P(cnt).Dp = P(cnt-1).Dp;
% P(cnt).vf = P(cnt-1).vf;
% v = [0.5 -1]; v = v/norm(v);
% p = ((X-n/2)*v(1) + (Y-n/2)*v(2));
% mask{cnt} = p < n*t & p > -n*t;
% dirs{cnt} = cat(3,-v(2)*ones(size(X)),v(1)*ones(size(X)));
% cnt = cnt + 1;


% % 
% % 
% P(3).Di = 2;
% P(3).Da = 1;
% P(3).Dp = 0.3;
% P(3).vf = 0.5;
% t = 0.2;
% v = [3 -3]; v = v/norm(v);
% p = ((X-n/2)*v(1) + (Y-n/2)*v(2));
% mask{3} = p < n*t & p > -n*t;
% dirs{3} = cat(3,-v(2)*ones(size(X)),v(1)*ones(size(X)));

signal = zeros([ size(X) size(ten,3)]);
for k = 1:length(mask),
    for j = 1:size(ten,3),
        E = ten(1,1,j)*dirs{k}(:,:,1).*dirs{k}(:,:,1) + 2*ten(1,2,j)*dirs{k}(:,:,1).*dirs{k}(:,:,2) + ten(2,2,j)*dirs{k}(:,:,2).*dirs{k}(:,:,2);
        Eorth = ten(1,1,j) + ten(2,2,j) + ten(3,3,j) - E;
        signal(:,:,j) = squeeze(signal(:,:,j)) + mask{k}.*(P(k).vsw + exp(-P(k).Di.*E).*P(k).vf + (1-P(k).vf-P(k).vsw).*exp(-P(k).Da.*E-P(k).Dp.*Eorth));
    end;
end;

fd = sum(cat(3,mask{:}),3);
signal = signal ./ (eps+repmat(fd,[ 1 1 size(ten,3)]));

themask = sum(cat(3,mask{:}),3)>0;

% 
% 
% figure(1);
% for k = 1:length(mask),
%     subplot(2,2,k);
%     imagesc(mask{k})
% end;
% figure(2);
% 
% imagesc(sum(cat(3,mask{:}),3));
szsig = size(signal);
signal = reshape(signal,[size(signal,1)*size(signal,2) size(signal,3)]);
Eiso = squeeze(ten(1,1,:) + ten(2,2,:) + ten(3,3,:));
signal(not(themask(:)),:) = repmat(exp(-Eiso*Dcsf)',[sum(not(themask(:))) 1]);
signal = reshape(signal,szsig);
%themask = themask*0 + 1;

signal = repmat(permute(signal,[1 2 4 3]),[1 1 d 1]);
signal = abs(signal + nz*(randn(size(signal)) +i*randn(size(signal))));

themask = repmat(themask,[1 1 d]);
%  
%  themask([1 end],:,:) = 0;
%  themask(:,[1 end],:) = 0;
  themask(:,:,[1 end]) = 0;
%  
% 
% 





function K = smker(ten)


%%
clear dir
for k = 1:size(ten,3),
    [U D] = eigs(ten(:,:,k));
    dir(:,k) = U(:,1);
end;

dir = [dir  -dir];

n = computeNeighbors(dir');











