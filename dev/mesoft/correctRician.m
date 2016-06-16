



x2dl = @(x,y,sigma,n) 1/sigma^2 * x.^2 - log(besseli(n/2-1,y.*x/sigma^2)) + (n/2-1)*log(y./x);


%%
x2dl = @(x,y,sigma,n) 1/sigma^2 * x.^2*(1+0*y') - log(besseli(n/2-1,x*y'/sigma^2)) + (n/2-1)*log((1./x)*y');


x = 0.01:0.0051:1; 
sigma = 0.05; 
[~,idx] = min(x2dl(x',x'*2,sigma,12)); 
y = x(idx);
[~,idx2] = sort(y); 
plot(x,y(idx)) 
axis([0 1 0 1])