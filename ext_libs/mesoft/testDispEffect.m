function testDispEffect
%%
ten = evalin('base','hr.user.bTensor')/1000;
for k = 1:size(ten,3),
    [U D] = eigs(ten(:,:,k));
    scheme(:,k) = sqrt(D(1,1))*U(:,1);
end;

D = 2; 
tries = 100;
n = randn(3,tries);
n = n ./ repmat(sqrt(sum(n.^2)),[3 1]);
a = 1-rand(tries,1)*0.2;

for k = 1:tries,
    S(:,k) = dispstick(D, a(k), n(:,k) ,scheme, 'gaussian');
end;
b = squeeze(round((ds.hr.user.bTensor(1,1,:)+ds.hr.user.bTensor(2,2,:)+ds.hr.user.bTensor(3,3,:))));

buni = unique(round(b/100))*100;
  
    for k = 2:length(buni),
            bidx = b == buni(k);                  
            [M shidx] = SH(dirs(:,bidx),lmax);
            pM = M*S(bidx); 
            pp = M*Smeas(bidx); 
            for l = 2:length(shidx)
                pMP(l-1,k-1) = pp(shidx{l})'*pM(shidx{l});
                pPred2(l-1,k-1) = sum(pM(shidx{l}).^2);
                pMeas2(l-1,k-1) = sum(pp(shidx{l}).^2);
            end;
        end;

    
function [M idx] = SH(n,lmax)
n = n ./ repmat(sqrt(sum(n.^2)),[3 1]);

M = [];
for l = 0:2:lmax,
    m = legendre(l,n(3,:),'sch');
    n2 = n(1,:)+i*n(2,:);
    n2 = n2 ./ abs(n2);
    m = m.* (repmat(n2,[size(m,1) 1]).^repmat((0:l)',[1 size(n2,2)]))*sqrt(2*l+1);
    idx1 = size(M,1);
    M = [M ; m(1,:) ; real(m(2:end,:)) ; imag(m(2:end,:))];
    idx2 = size(M,1);
    idx{l/2+1} = idx1+1:idx2;
end;

M = M/sqrt(size(n,2));

