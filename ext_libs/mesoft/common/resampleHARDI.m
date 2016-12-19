function mr = resampleHARDI(mr,bDir)


   bTensor = mr.user.bTensor;
  
   bval = sum(mr.user.bDir.^2);
   b0idx = bval == 0;
   didx = bval>0;

   bvalout = sum(bDir.^2);
   b0idxout = bvalout == 0;
   didxout = bvalout>0;
   
   
   sz = size(mr.dataAy);
   
   
   b0img = mean(mr.dataAy(:,:,:,b0idx),4);
      
   sinterp = sphereInterpolLUT(mr.user.bDir(:,didx)');
   
   
   sig = single(permute(mr.dataAy(:,:,:,didx),[4 1 2 3]));   
   
   tdirfield = repmat(single(bDir(:,didxout)),[1 1 prod(sz(1:3))]);
   sig = evalSinterp(tdirfield(:,:,:),sig(:,:),sinterp);
   sig = reshape(sig,[sum(didxout) sz(1:3)]);
   

   numoutDir = size(bDir,2);
   mr.dataAy = zeros(sz(1),sz(2),sz(3),numoutDir,'single');
   mr.dataAy(:,:,:,didxout) = permute(sig,[2 3 4 1]);
   mr.dataAy(:,:,:,b0idxout) = repmat(b0img,[1 1 1 sum(b0idxout)]);
   mr.user.bDir = bDir;
   mr.user.bTensor = [];
   for k = 1:size(bDir,2),
       mr.user.bTensor(:,:,k) = bDir(:,k)*bDir(:,k)';
   end;
%    
% 
% function mr = resampleHARDI(mr,bDir)
% 
%    bTensor = mr.user.bTensor;
%   
%    bval = sum(mr.user.bDir.^2);
%    b0idx = bval == 0;
%    didx = bval>0;
% 
%    bvalout = sum(bDir.^2);
%    b0idxout = bvalout == 0;
%    didxout = bvalout>0;
%    
%    
%    sz = size(mr.dataAy);
%    
%    
%    b0img = mean(mr.dataAy(:,:,:,b0idx),4);
%    resampler = Resample(bTensor(:,:,didx),bDir(:,didxout),40,0.00000006);      
%    sig = permute(mr.dataAy(:,:,:,didx),[4 1 2 3]);   
%    sig = reshape(resampler*sig(:,:),[sum(didxout) sz(1:3)]);
%    
% 
%    numoutDir = size(bDir,2);
%    mr.dataAy = zeros(sz(1),sz(2),sz(3),numoutDir);
%    mr.dataAy(:,:,:,didxout) = permute(sig,[2 3 4 1]);
%    mr.dataAy(:,:,:,b0idxout) = repmat(b0img,[1 1 1 sum(b0idxout)]);
%    mr.user.bDir = bDir;
%    mr.user.bTensor = [];
%    for k = 1:size(bDir,2),
%        mr.user.bTensor(:,:,k) = bDir(:,k)*bDir(:,k)';
%    end;
%    
       
function [odf odfmean] = Resample(pos,posout,L,alpha)

N = size(pos,3);

SH = cpmSH(pos,L);
SH = vertcat(SH{:});
SH = SH ./ sqrt((ones(size(SH,1),1)*sum((abs(SH)).^2)));
SHout = cpmSH(posout,L);
SHout = vertcat(SHout{:});
SHout = SHout ./ sqrt((ones(size(SHout,1),1)*sum((abs(SHout)).^2)));

for k = 1:2:L
    regu{k} = ones(2*k-1,1)*(k-1)*k;
end;

regu = vertcat(regu{:});

odfsh = inv(SH*SH'+alpha*N*diag(regu))*SH;

FRT = SHout'*odfsh;

odf = real(FRT);
odfmean = real(odfsh(1,:,:,:));

 



