
function [samples] = mf_sample(mu)


kappa = eps+sqrt(sum(mu.^2));
mu = mu./ repmat(kappa,[size(mu,1) 1]);

kappa = kappa(:)';
kappa(kappa>200) = 200;

d = size(mu,1);
n = size(mu,2);


iso = kappa < (d-1)/2;
kappa(iso) = (d-1)/2;


b = (d-1)./( 2*kappa + sqrt((4*kappa.^2 - (d-1)^2))); 
x0 = (1 - b)./(1 + b); 
c = kappa.*x0 + (d - 1).*log(1-x0.^2);

m = size(b,2);
idx = 1:m;
W = zeros(1,m);

cnt = 1;
while m>0  ,
    cnt = cnt + 1;
    Z = betarnd( (d-1)/2, (d-1)/2, 1,m);
    Wtmp = (1 - (1 + b).*Z) ./ (1 - (1 - b).*Z+eps);
    valid = kappa.*Wtmp + (d - 1).*log(1 - x0.*Wtmp) - c > log(rand(1,m));
    W(idx(valid)) = Wtmp(valid);
    c = c(not(valid));
    b = b(not(valid));
    x0 = x0(not(valid));
    kappa = kappa(not(valid));
    idx = idx(not(valid));
    m = length(idx);
    if cnt > 100,
        break;
    end;
end;

V = randn(d,n);
Vorth = ( V - mu.*repmat(sum(V.*mu),[d 1]) );
Vorth = Vorth ./ repmat(sqrt(sum(Vorth.^2)),[d 1]);
samples = repmat(sqrt(abs(1-W.^2)),[d 1]) .* Vorth + mu.*repmat(W,[d 1]);



