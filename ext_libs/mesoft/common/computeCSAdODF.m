function mr = computeCSAdODF(mr,numdirs,Lmax,penalty)


   dwidirs = load('dwidirections');
   outdir = eval(sprintf('dwidirs.dirs%i',numdirs));

   bTensor = mr.user.bTensor;
  
   bval = sum(mr.user.bDir.^2);
   b0idx = bval == 0;
   didx = bval>0;

  
   
   sz = size(mr.dataAy);
   
   
   b0img = mean(mr.dataAy(:,:,:,b0idx),4);
   
   signal = mr.dataAy(:,:,:,didx) ./ repmat(b0img,[1 1 1 sum(didx)]);
   signal = permute(signal,[4 1 2 3]);   

   ss = 0.999;       
   signal = (signal>ss)*ss + (signal<=ss).*signal;
   signal = log(-log(signal));
   [frt frtmean] = FunkRadonTransLaplace(bTensor(:,:,didx),outdir,Lmax,penalty);      
   sz = size(signal);
   sz(1) = size(frt,1);
   signal = 1+single(reshape(frt*signal(:,:),sz));
   
   
   mr.dataAy = permute(signal,[2 3 4 1]); clear signal;
   mr.user.b0avg = b0img;
   mr.user.mFOD = mean(mr.dataAy,4);
   mr.user.mask = (mr.user.b0avg*0+1)>0;
   mr.user.sym = true;
   mr.user.parameters.Lmax = Lmax;
   mr.user.parameters.penalty = penalty;
   mr.user.parameters.algorithm = 'CSAAganj';
   mr.user.bDir = outdir;
   
   
   


       
function [odf odfmean] = Resample(pos,posout,L,alpha)

L = L +1;
N = size(pos,3);

SH = cpmSH(pos,L);
SH = vertcat(SH{:});
SH = SH ./ sqrt((ones(size(SH,1),1)*sum((abs(SH)).^2)));
SHout = cpmSH(posout,L);
SHout = vertcat(SHout{:});
SHout = SHout ./ sqrt((ones(size(SHout,1),1)*sum((abs(SHout)).^2)));

for k = 1:2:L
    regu{k} = ones(2*k-1,1)*((k-1)*k)^2;
end;

regu = vertcat(regu{:});

odfsh = inv(SH*SH'+alpha*N*diag(regu))*SH;

FRT = SHout'*odfsh;

odf = real(FRT);
odfmean = real(odfsh(1,:,:,:));

 



