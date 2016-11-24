    
function [est cmap] = compMesoParams(hr);




%%
keyboard
sz = size(hr.dataAy);
ten = hr.user.bTensor/1000;
for k = 1:size(ten,3),
    
    [U D] = eigs(ten(:,:,k));
    [~,ix]=sort(D(logical(eye(length(D)))));
    scheme(:,k) = sqrt(D(1,1))*U(:,ix(3));
    
end;



sigma = 1.75;
ng = ceil(sigma*3); ng = ng + mod(ng,2) + 1;
gau = fspecial('gaussian',[ng 1],sigma); [gx gy gz] = ndgrid(gau); gau = gx.*gy.*gz;
sz = size(hr.dataAy); sz=sz(1:3);
gau = padarray(gau,floor((sz-ng)/2),'post');
gau = padarray(gau,floor((sz-ng)/2)-mod(sz,2) +1,'pre');
gau = fftn(ifftshift(gau));
sm = @(x) real(ifftn(fftn(x).*gau));   


b = sum(scheme.^2);
buni = unique(round(b*10))/10;
b0idx = round(b*10)/10==0;
b0 = mean(hr.dataAy(:,:,:,b0idx),4);
err = std(hr.dataAy(:,:,:,b0idx),[],4);

valid = hr.dataAy(:,:,:,1)>0; valid = valid(:);
%SNR = b0/ mean(err(valid));
SNR = b0./ (eps+err) ;

%  
%  for k = 1:size(hr.dataAy,4);
%      hr.dataAy(:,:,:,k) = sm(hr.dataAy(:,:,:,k));
%  end;
 
SNR = sm(SNR);




modelstruc = evalin('base','modelstruc');
S = reshape(hr.dataAy,[size(hr.dataAy,1)*size(hr.dataAy,2)*size(hr.dataAy,3) size(hr.dataAy,4)]);
S = S./(eps+repmat(b0(:),[1 size(hr.dataAy,4)]));  

[est cmap] = modelstruc.apply(SNR(:),S,@(x)(x));
est = reshape(est,[sz(1:3) 6]);
cmap = reshape(cmap,sz(1:3));

est(est>5) = 3;
est(est<0) = 0;


est(:,:,:,end+1) = SNR;

mask = (b0>150) & (cmap>0)  ;%.*(est(:,:,:,4)>0.05));
%%
mask = bwareaopen(mask,20,6);
mask = not(bwareaopen(not(mask),10,6));

%%
est = est.*repmat(mask,[1 1 1 size(est,4)]);



% est(isnan(est(:))) = 0;
% for k = 1:size(est,4);
%     est(:,:,:,k) = sm(est(:,:,:,k));
% end;

return;







    
    
    
function [M idx] = SH(n,lmax)
n = n ./ repmat(sqrt(sum(n.^2)),[3 1]);

M = [];
for l = 0:2:lmax,
    m = legendre(l,n(3,:),'sch');
    n2 = n(1,:)+i*n(2,:);
    n2 = n2 ./ abs(n2);
    m = m.* (repmat(n2,[size(m,1) 1]).^repmat((0:l)',[1 size(n2,2)]));
    idx1 = size(M,1);
    M = [M ; m(1,:) ; real(m(2:end,:)) ; imag(m(2:end,:))];
    idx2 = size(M,1);
    idx{l/2+1} = idx1+1:idx2;
end;

M = M/sqrt(size(n,2));

            





