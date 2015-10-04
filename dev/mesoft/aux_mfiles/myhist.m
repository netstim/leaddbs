
function [h b] = myhist(d,w,ran,n)
d = double(floor((d-ran(1))/(ran(2)-ran(1))*(n-1) +1));
w = w(d>=1 & d<= n);
d = d(d>=1 & d<= n);
h = sparse([d(:);n],ones(length(d(:))+1,1),double([w(:) ;eps]));
h = full(h);
b = (1:n)'/n*(ran(2)-ran(1)) + ran(1);
