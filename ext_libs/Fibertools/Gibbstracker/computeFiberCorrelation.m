function [coeff meanval_sq] = computeFiberCorrelation(bDir,bD)
n = size(bDir,2); 

C = abs(bDir'*bDir);

Q = exp(-bD * C.^2);
Q = Q - repmat(mean(Q),[size(Q,1) 1 ]);


P = Q*Q;
[alpha idx] = sort(C(:));
beta = P(idx);

nfac = sqrt(beta(end));
beta = beta / nfac^2;
Q = Q /nfac;
meanval_sq = 0; %sum(Q(1,:))^2/n;


%beta = beta - meanval_sq;% +chempot_2nd;
%beta = beta - min(Q(:));

T = [alpha.^0 alpha.^2 alpha.^4 alpha.^6];
coeff = pinv(T)*beta;

x = (0:0.01:1)';
approx = [x.^0  x.^2 x.^4 x.^6]*coeff;

%figure(11); plot(acos(alpha),beta,acos(x),[x.^0  x.^2 x.^4 x.^6]*coeff);

