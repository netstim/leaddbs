%% mises

e = '-(a.^(3./2).*(6.*exp(-1./a) - 6.*exp(1./a) + (6.*exp(1./a))./a - (2.*exp(1./a))./a.^2 + (6.*exp(-1./a))./a + (2.*exp(-1./a))./a.^2))./(4.*(1./a).^(1./2).*(exp(1./a)./2 - exp(-1./a)./2))';


a = 0:0.01:2;
f = eval(e);
f(1) = 1;
plot(a,f)
kappa = 0.001;
n = 8; 
M =( repmat(f,[n 1]).^repmat((0:(n-1))',[1 size(f,2)]) )';

as = exp(M*pinv(M)*log(a'+kappa))-kappa
 plot(as,f,'r',a,f,'b')
 
 %%
%alpha =  (pinv(M)*log(a'+kappa))'
%alpha = [0.0011   -0.0306    0.2841   -1.5668    4.7801   -7.9515    6.6931   -2.2164]*10^3;
alpha = [ 0.0011   -0.0290    0.2473   -1.2354    3.3866   -5.0763    3.8899   -1.1908]*10^3;
SH2ic = @(f) 1./(exp((repmat(f,[n 1]).^repmat((0:(n-1))',[1 size(f,2)]) )'*alpha')-kappa);
t = 0.01:0.01:0.98; plot(t,SH2ic(t))

%% watson

e = '(((3.*a.*exp(1./a))./2 - (pi.^(1./2).*erfi((1./a).^(1./2)).*((3.*a)./4 + 1./2))./(1./a).^(1./2)).*(1./a).^(1./2))./(pi.^(1./2).*erfi((1./a).^(1./2)))';


a = 0.01:0.01:6;
f = eval(e);
f(1) = 1;
plot(a,f)
kappa = 0.001;
n = 8; 
M =( repmat(f,[n 1]).^repmat((0:(n-1))',[1 size(f,2)]) )';

as = exp(M*pinv(M)*log(a'+kappa))-kappa
plot(as,f,'r',a,f,'b')

alpha =  (pinv(M)*log(a'+kappa))'
%alpha = [0.0021   -0.0298    0.1735   -0.6539    1.4737   -1.9221    1.3324   -0.3803]*10^3;
SH2c = @(f) 1./(exp((repmat(f,[n 1]).^repmat((0:(n-1))',[1 size(f,2)]) )'*alpha')-kappa);
t = [0.0:0.01:1];
plot(t,SH2c(t))

