function [odf odfmean] = FunkRadonTrans(pos,posout,L,alpha)

N = size(pos,3);
SH = cpmSH(pos,L);
SH = vertcat(SH{:});
SHout = cpmSH(posout,L);
SHout = vertcat(SHout{:});

for k = 1:2:L
    regu{k} = ones(2*k-1,1)*(k-1)*k;
    Leg = legendre(k-1,0,'norm');
    funkradon{k} =ones(2*k-1,1)*Leg(1)*3^(k-1)/factorial(k-1);         
end;

funkradon{1} = 0;

regu = vertcat(regu{:});
funkradon = vertcat(funkradon{:});

odfsh = inv(SH*SH'+alpha*N*diag(regu))*SH;

FRT = SHout'*diag(funkradon)*odfsh;

odf = real(FRT);
odfmean = real(odfsh(1,:,:,:));



return
% 
% U = [ -1/sqrt(2) 0 0 ; ...
%        0   -1/sqrt(2)  0 ; ...  
%        0    0      1 ; ... 
%        0    0      0 ];
% udir = single(U*dir);
% SH{1} = single([ones(1,N) ; zeros(1,N)]);
% SH{2} = single(udir);
% for k = 3:L,
%     SH{k} = STmultiply(udir,SH{k-1},k-1,1/ClebschGordan(1,k-2,k-1,0,0,0));
% end;
% for k = 1:L,
%     SH{k} = sqrt(2*(k-1)+1) * SH{k}(1:(end-1),:);
% end;
    

   








