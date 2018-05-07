function S = dispstick(D, a, n,sch,type)

if length(size(sch)) == 3,
    for k = 1:size(sch,3),
        [U S] = eigs(sch(:,:,k));
        scheme(:,k) = U(:,1) *sqrt(S(1,1));
    end;
elseif size(sch,1) == 2,
    n = [1 0 0]';
    scheme = [sqrt(sch(1,:)) ; sqrt(sch(2,:)) ; sch(2,:)*0];
else
    scheme = sch;
end;


if strcmp(type,'gaussian')==1,
    avg = 500;
    a = -log(a)/2;
    nz = repmat(n,[1 avg]) + a * randn(3,avg);
    nz = nz ./ repmat(sqrt(sum(nz.^2)),[3 1]);
    qxx = (scheme'*nz).^2;    
    S = mean(exp(-D*qxx),2);
elseif strcmp(type,'mises')==1,
    avg = 500;
        
    kappa = 0.001;
    alpha = [ 0.0011   -0.0290    0.2473   -1.2354    3.3866   -5.0763    3.8899   -1.1908]*10^3;
    SH2ic = @(f) 1./(exp((repmat(f,[8 1]).^repmat((0:(8-1))',[1 size(f,2)]) )'*alpha')-kappa);    
    cp = SH2ic(a);
    
    nz = mf_sample(repmat(n,[1 avg])*cp);
    qxx = (scheme'*nz).^2;
    S = mean(exp(-D*qxx),2);
elseif strcmp(type,'poisson')==1,
    a = sqrt(a);
    S = poissonDisp((scheme'*n).^2,sum(scheme.^2),D,a);    
elseif strcmp(type,'heat')==1,
    a = -log(a)/6;
    S = heatDisp((scheme'*n).^2,sum(scheme.^2),D,a);    
elseif strcmp(type,'watson')==1,
    S = watsonDisp((scheme'*n).^2,sum(scheme.^2),D,a);    
end;



function S = poissonDisp(qxx,b,D,lam)
lmax = 50;
L = (0:lmax);
f = lam.^L ;
S = SHDisp(qxx,b,D,f);


function S = heatDisp(qxx,b,D,lam)
lmax = 50;
L = (0:lmax);
f = exp(-lam*L.*(L+1));
S = SHDisp(qxx,b,D,f);


function S = watsonDisp(qxx,b,D,lam)

watsonlut = evalin('base','watsonlut');

k = 0.001;
alpha = [0.0026   -0.0449    0.3326   -1.4291    3.4388   -4.5957    3.1820   -0.8907]*10^3;
SH2c = @(f) 1./(exp((repmat(f,[8 1]).^repmat((0:(7))',[1 size(f,2)]) )'*alpha')-k);
kappa = SH2c(lam);
kappa(kappa<0) = 0;
lmax = 14;
%%
f(1) = 1;
for l = 2:2:lmax,
    f(l+1) = interp1(watsonlut.k,watsonlut.lut(l/2,:),kappa,'linear','extrap') ;
end;
S = SHDisp(qxx,b,D,f);







function S = SHDisp(qxx,b,D,f)

buni = unique(round(b*10));
t = -1:0.001:1;
lmax = length(f)-1;
L = (0:lmax);
p = myleg(lmax,t') .* repmat(sqrt((2*L+1)),[size(t,2) 1]) /sqrt(length(t));
S = ones(length(b),1);
for k = 2:length(buni),
    bD = D*buni(k)/10;
    idx = round(b*10)==buni(k);
    P = ((exp(-bD*t.^2)*p).*f)*p';
    it = sqrt(qxx(idx)/buni(k)*10); it(it>1) = 1; it(it<0) = 0;
    S(idx) = interp1(t,P,it);    
end;




return;
    






function p = myleg(n,x);
if n == 0,
    p = x*0+1;
    return;
end
if n == 1
    p = [x*0+1 x];
    return;
end;
p = zeros(size(x,1),n+1);
p(:,1:2) = [x*0+1 x];
for k = 2:n,
    p(:,k+1) = ((2*k-1)*x.*p(:,k) - (k-1)*p(:,k-1))/k;    
end;


function [C] = WatsonSHCoeff(k)
% function [C, D] = WatsonSHCoeff(k)
% Computes the spherical harmonic (SH) coefficients of the Watson's
% distribution with the concentration parameter k (kappa) up to the 12th order
% and the derivatives if requested.
%
% Truncating at the 12th order gives good approximation for kappa up to 64.
%
% INPUTS:
%
% k should be an array of positive numbers, specifying a set of
% concentration parameters for the Watson's distribution.
%
% OUTPUTS:
%
% C will be a 2-D array and each row contains the SH coefficients of the
% orders 0, 2, 4, ..., to 2n for the parameter in the corresponding row in
% k.
%
% Note that the SH coefficients of the odd orders are always zero.
%
% D will be the 1st order derivative of C.
%
% author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%

large = find(k>30);
exact = find(k>0.1);
approx = find(k<=0.1);
% Necessary to make matlab happy when k is a single value
exact = exact(:);
approx = approx(:);
large = large(:);

% The maximum order of SH coefficients (2n)
n = 6;

% Computing the SH coefficients
C = zeros(length(k),n+1);

% 0th order is a constant
C(:,1) = 2*sqrt(pi);

% Precompute the special function values
sk = sqrt(k(exact));
sk2 = sk.*k(exact);
sk3 = sk2.*k(exact);
sk4 = sk3.*k(exact);
sk5 = sk4.*k(exact);
sk6 = sk5.*k(exact);
sk7 = sk6.*k(exact);
k2 = k.^2;
k3 = k2.*k;
k4 = k3.*k;
k5 = k4.*k;
k6 = k5.*k;
k7 = k6.*k;

erfik = erfi(sk);
ierfik = 1./erfik;
ek = exp(k(exact));
dawsonk = 0.5*sqrt(pi)*erfik./ek;

% for large enough kappa
C(exact,2) = 3*sk - (3 + 2*k(exact)).*dawsonk;
C(exact,2) = sqrt(5)*C(exact,2).*ek;
C(exact,2) = C(exact,2).*ierfik./k(exact);

C(exact,3) = (105 + 60*k(exact) + 12*k2(exact)).*dawsonk;
C(exact,3) = C(exact,3) -105*sk + 10*sk2;
C(exact,3) = .375*C(exact,3).*ek./k2(exact);
C(exact,3) = C(exact,3).*ierfik;

C(exact,4) = -3465 - 1890*k(exact) - 420*k2(exact) - 40*k3(exact);
C(exact,4) = C(exact,4).*dawsonk;
C(exact,4) = C(exact,4) + 3465*sk - 420*sk2 + 84*sk3;
C(exact,4) = C(exact,4)*sqrt(13*pi)/64./k3(exact);
C(exact,4) = C(exact,4)./dawsonk;

C(exact,5) = 675675 + 360360*k(exact) + 83160*k2(exact) + 10080*k3(exact) + 560*k4(exact);
C(exact,5) = C(exact,5).*dawsonk;
C(exact,5) = C(exact,5) - 675675*sk + 90090*sk2 - 23100*sk3 + 744*sk4;
C(exact,5) = sqrt(17)*C(exact,5).*ek;
C(exact,5) = C(exact,5)/512./k4(exact);
C(exact,5) = C(exact,5).*ierfik;

C(exact,6) = -43648605 - 22972950*k(exact) - 5405400*k2(exact) - 720720*k3(exact) - 55440*k4(exact) - 2016*k5(exact);
C(exact,6) = C(exact,6).*dawsonk;
C(exact,6) = C(exact,6) + 43648605*sk - 6126120*sk2 + 1729728*sk3 - 82368*sk4 + 5104*sk5;
C(exact,6) = sqrt(21*pi)*C(exact,6)/4096./k5(exact);
C(exact,6) = C(exact,6)./dawsonk;

C(exact,7) = 7027425405 + 3666482820*k(exact) + 872972100*k2(exact) + 122522400*k3(exact)  + 10810800*k4(exact) + 576576*k5(exact) + 14784*k6(exact);
C(exact,7) = C(exact,7).*dawsonk;
C(exact,7) = C(exact,7) - 7027425405*sk + 1018467450*sk2 - 302630328*sk3 + 17153136*sk4 - 1553552*sk5 + 25376*sk6;
C(exact,7) = 5*C(exact,7).*ek;
C(exact,7) = C(exact,7)/16384./k6(exact);
C(exact,7) = C(exact,7).*ierfik;

% for very large kappa
if size(large,1) > 0
  lnkd = log(k(large)) - log(30);
  lnkd2 = lnkd.*lnkd;
  lnkd3 = lnkd2.*lnkd;
  lnkd4 = lnkd3.*lnkd;
  lnkd5 = lnkd4.*lnkd;
  lnkd6 = lnkd5.*lnkd;
  C(large,2) = 7.52308 + 0.411538*lnkd - 0.214588*lnkd2 + 0.0784091*lnkd3 - 0.023981*lnkd4 + 0.00731537*lnkd5 - 0.0026467*lnkd6;
  C(large,3) = 8.93718 + 1.62147*lnkd - 0.733421*lnkd2 + 0.191568*lnkd3 - 0.0202906*lnkd4 - 0.00779095*lnkd5 + 0.00574847*lnkd6;
  C(large,4) = 8.87905 + 3.35689*lnkd - 1.15935*lnkd2 + 0.0673053*lnkd3 + 0.121857*lnkd4 - 0.066642*lnkd5 + 0.0180215*lnkd6;
  C(large,5) = 7.84352 + 5.03178*lnkd - 1.0193*lnkd2 - 0.426362*lnkd3 + 0.328816*lnkd4 - 0.0688176*lnkd5 - 0.0229398*lnkd6;
  C(large,6) = 6.30113 + 6.09914*lnkd - 0.16088*lnkd2 - 1.05578*lnkd3 + 0.338069*lnkd4 + 0.0937157*lnkd5 - 0.106935*lnkd6;
  C(large,7) = 4.65678 + 6.30069*lnkd + 1.13754*lnkd2 - 1.38393*lnkd3 - 0.0134758*lnkd4 + 0.331686*lnkd5 - 0.105954*lnkd6;
end

% for small kappa
C(approx,2) = 4/3*k(approx) + 8/63*k2(approx);
C(approx,2) = C(approx,2)*sqrt(pi/5);

C(approx,3) = 8/21*k2(approx) + 32/693*k3(approx);
C(approx,3) = C(approx,3)*(sqrt(pi)*0.2);

C(approx,4) = 16/693*k3(approx) + 32/10395*k4(approx);
C(approx,4) = C(approx,4)*sqrt(pi/13);

C(approx,5) = 32/19305*k4(approx);
C(approx,5) = C(approx,5)*sqrt(pi/17);

C(approx,6) = 64*sqrt(pi/21)*k5(approx)/692835;

C(approx,7) = 128*sqrt(pi)*k6(approx)/152108775;


