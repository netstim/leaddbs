function pf = fdist(f,v1,v2)
%FDIST( F, v1, v2) returns Q(F|v1,v2), the probability
%   of observing a value of F or greater in an
%   F-distribution with v1 and v2 degrees of freedom.
%
%   e.g., in Abramowitz & Stegun's Table 26.9, the value
%   F for which Q(F|v1,v2)=0.01 when v1=3, v2=20 is 4.94;
%   this can be verified by confirming that
%   FDIST( 4.94, 3, 20 ) returns a probability of 0.01, or 1%.
%
%   To find the value 4.94, use FDISTINV.
%   To apply the F-test, use FTEST_MODELS
%   also see FTEST_BACKGROUND and FTEST_EXAMPLE

% References: Press et al., Numerical Recipes,
%   Cambridge, 1986;
% Abramowitz & Stegun, Handbook of Mathematical
%   Functions, Dover, 1972.

% Peter R. Shaw, Woods Hole Oceanographic Institution
% Woods Hole, MA 02543
% (508) 457-2000 ext. 2473  pshaw@aqua.whoi.edu
% March, 1990

% ^ Calls functions BETAI, BETACF and GAMMLN  ^

a = v2 ./ 2;
b = v1 ./ 2;
x = v2 ./ ( v2 + v1 .* f) ;
if(a<0),
 pf = 1.0;
else
 pf = local_betai(a,b,x);
end

function bi= local_betai(a,b,X)
%BETAI  Incomplete Beta function.
%  BETAI(a,b,x) returns the Incomplete Beta function Ix(a,b)
%  for every element of x.  Parameters a and b must be scalars.

% Peter R. Shaw, Woods Hole Oceanographic Institution
% Woods Hole, MA 02543
% (508) 457-2000 ext. 2473  pshaw@aqua.whoi.edu

% Converted from the Fortran subroutine "BETAI" in:
% Numerical Recipes, Press et al., Cambridge, 1986.

% ^ Calls functions BETACF and GAMMLN  ^

[m,n]=size(X);
bi=zeros(m,n);
%sorry for this non-vectorized loop, folks:
for i=1:m,
  for j=1:n,
    x=X(i,j);
    if x<0 || x>1,
      error('bad argument x in BETAI')
    end
    if x==0.0 || x==1.0,
      bt=0.0;  % Factor in front of continued fraction
    else
      bt=exp(local_gammln(a+b)-local_gammln(a)-local_gammln(b) ...
           +a*log(x)+b*log(1.0-x));
    end
    if x<(a+1)/(a+b+2), %use continued fraction directly.
      bi(i,j) = bt*local_betacf(a,b,x)/a;
    else
%     Use continued fraction after making
%     symmetry transformation:
      bi(i,j) = 1.0-bt*local_betacf(b,a,1.0-x)/b;
    end
  end
end

function bcf= local_betacf(a,b,x)
%local_betacf(a,b,x) is a continued fraction representation
%  of the elements of x;
%  used in evaluating the incomplete Beta function BETAI.

% Peter R. Shaw, Woods Hole Oceanographic Institution
% Woods Hole, MA 02543
% (508) 457-2000 ext. 2473  pshaw@aqua.whoi.edu

% Converted from the Fortran subroutine "local_betacf" in:
% Numerical Recipes, Press et al., Cambridge, 1986.

% ^ Calls no other routines ^

[mmx,nnx]=size(x);
bcf=zeros(mmx,nnx);
itmax=100; epsilon=3.e-7;

%sorry for this non-vectorized loop, folks:
for i=1:mmx,
  for j=1:nnx,
    am=1; bm=1; az=1;
    qab=a+b;
    qap=a+1;
    qam=a-1;
    bz=1-qab*x(i,j)/qap;
    aold=0; az=1; m=0;
    while abs(az-aold) > epsilon*abs(az),
       m=m+1;
       if m>itmax,
          error('(local_betacf): a or b too big or itmax too small');
       end
       em=m;
       tem=em+em;
       d=em*(b-m)*x(i,j)/((qam+tem)*(a+tem));
       ap=az+d*am;
       bp=bz+d*bm;
       d=-(a+em)*(qab+em)*x(i,j)/((a+tem)*(qap+tem));
       app=ap+d*az;
       bpp=bp+d*bz;
       aold=az;
       am=ap/bpp;
       bm=bp/bpp;
       az=app/bpp;
       bz=1.0 ;
    end
    bcf(i,j)=az;
  end
end


function gl=local_gammln(xx)
%GAMMLN   Natural log of the complete Gamma function.
%   GAMMLN(X) returns the log of the gamma function
%   for every element of X.

%   Useful in formulas involving, e.g., gamma(x)/gamma(y)
%     for large x and y.

% Peter R. Shaw, Woods Hole Oceanographic Institution
% Woods Hole, MA 02543
% (508) 457-2000 ext. 2473  pshaw@aqua.whoi.edu

% Converted from the Fortran subroutine "GAMMLN" in:
% Numerical Recipes, Press et al., Cambridge, 1986.

% ^ Calls no other routines ^

cof=[76.18009173, -86.50532033, 24.01409822, ...
-1.231739516, 0.120858003e-2, -0.536382e-5];
stp=2.50662827465;
[mxx,nxx]=size(xx);
one=ones(mxx,nxx);
half=0.5*one; fpf=5.5*one;
x=xx-one;
tmp=x+fpf;
tmp=(x+half).*log(tmp)-tmp;
ser=one;
for j=1:6,
  x=x+one;
  ser=ser+cof(j)./x;
end
gl=tmp+log(stp*ser);
