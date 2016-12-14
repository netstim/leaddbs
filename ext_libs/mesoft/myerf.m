function res = myerf(x)

res = x;
idx = abs(imag(x))>eps;
res(idx) = erfi(imag(x(idx)))*1.i;
res(not(idx)) = erf(x(not(idx)));