
function hrsm = smoothAnisoRice(hr,axsig,rasig,ksz,sigma);


if sigma == 0,
    xx = sort(hr.dataAy(:));
    sigma = xx(round(length(hr.dataAy(:))*0.2));
end;


hr.dataAy = hr.dataAy/sigma;
hrsm = hr;
hrsm.dataAy = hr.dataAy*0;

bDir = hr.user.bDir;
sz = size(hr.dataAy);
[X Y Z] = ndgrid(-ksz:ksz);
R2 = X.^2 + Y.^2 + Z.^2;
isog = (exp(-R2/(2*axsig^2)));       

for k = 1:size(hr.dataAy,4)
   fprintf('.');
   fgauss = exp(- (R2-(X*bDir(1,k) +  Y*bDir(2,k) + Z*bDir(3,k)).^2) * (1/(2*rasig.^2) - 1/(2*axsig.^2))) .* isog;
   fgauss = fgauss/sum(fgauss(:));
   %fgauss = fgauss(:);
%    for z = 1+ksz:size(hr.dataAy,3)-ksz,
%        fprintf(',');
%        for y = 1+ksz:size(hr.dataAy,2)-ksz,
%             for x = 1+ksz:size(hr.dataAy,1)-ksz,
%                 a = hr.dataAy(x-ksz:x+ksz,y-ksz:y+ksz,z-ksz:z+ksz,k); 
%                 hrsm.dataAy(x,y,z,k) = RicianMean(a(:),fgauss);
%             end
%        end
%    end
    hrsm.dataAy(:,:,:,k) = smoothRice(double(hr.dataAy(:,:,:,k)),double(fgauss));


   fprintf('\n');       
end
hrsm.dataAy = single(hrsm.dataAy*sigma);

function X = RicianMean(a,w);

X = sum(a.*w)/sum(w(:));

for it = 1:100,
    R = besseli(1,X*a)./besseli(0,X*a);
    delta = sum(  w.*(X-a.*R) ) /sum( w.*((a.*R).^2 - (a.^2 - a.*R/X) + 1)) ;    
    X = X - delta;
    if delta < 10^-3,
        break;
    end;
end;


    
    
    

