function R = chiLogLik(S,M,n,sig)
% S - model signal
% M - measurement

S = double(S/sig);
M = double(M/sig);
lb = S.*M;
lb(lb<500) = log(besseli(n/2-1,lb(lb<500)));
R = sum ((S.^2 + M.^2)/2 - n/2*log(M) - (1-n/2)*log(S) -lb);