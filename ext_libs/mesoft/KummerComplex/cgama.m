function [gr,gi]=cgama(x,y,kf)

%  [RealPartResult,ImaginaryPartResult]=cgama(RealPart,ImaginaryPart,1)
%

%This program is a direct conversion of the corresponding Fortran program in
%S. Zhang & J. Jin "Computation of Special Functions" (Wiley, 1996).
%online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html
%
%Converted by f2matlab open source project:
%online: https://sourceforge.net/projects/f2matlab/
% written by Ben Barrowes (barrowes@alum.mit.edu)
%
%     ==========================================================
%     Purpose: This program computes the gamma function a(z)
%     or ln[a(z)] for a complex argument using
%     subroutine CGAMA
%     Input :  x  --- Real part of z
%     y  --- Imaginary part of z
%     KF --- Function code
%     KF=0 for ln[a(z)]
%     KF=1 for a(z)
%     Output:  GR --- Real part of ln[a(z)] or a(z)
%     GI --- Imaginary part of ln[a(z)] or a(z)
%     Examples:
%     x         y           Re[a(z)]           Im[a(z)]
%     --------------------------------------------------------
%     2.50      5.00     .2267360319D-01    -.1172284404D-01
%     5.00     10.00     .1327696517D-01     .3639011746D-02
%     2.50     -5.00     .2267360319D-01     .1172284404D-01
%     5.00    -10.00     .1327696517D-01    -.3639011746D-02
%     x         y          Re[lna(z)]         Im[lna(z)]
%     ---------------------------------------------------------
%     2.50      5.00    -.3668103262D+01     .5806009801D+01
%     5.00     10.00    -.4285507444D+01     .1911707090D+02
%     2.50     -5.00    -.3668103262D+01    -.5806009801D+01
%     5.00    -10.00    -.4285507444D+01    -.1911707090D+02
%     ==========================================================



%     for a complex argument
%     Input :  x  --- Real part of z
%     y  --- Imaginary part of z
%     KF --- Function code
%     KF=0 for ln[a(z)]
%     KF=1 for a(z)
%     Output:  GR --- Real part of ln[a(z)] or a(z)
%     GI --- Imaginary part of ln[a(z)] or a(z)
%     ========================================================
%
%
%
%

a=zeros(10,1);
x1=0.0;
pi=3.141592653589793d0;
a=[8.333333333333333d-02,-2.777777777777778d-03,7.936507936507937d-04,-5.952380952380952d-04,8.417508417508418d-04,-1.917526917526918d-03,6.410256410256410d-03,-2.955065359477124d-02,1.796443723688307d-01,-1.39243221690590d+00];

if (y == 0.0d0&&x == fix(x)&&x <= 0.0d0) ;
	gr=NaN;
	gi=NaN;
	return;
elseif (x < 0.0d0);
	x1=x;
	y1=y;
	x=-x;
	y=-y;
end;

x0=x;

if (x <= 7.0) ;
	na=fix(7-x);
	x0=x+na;
end;

z1=sqrt(x0.*x0+y.*y);
th=atan(y./x0);
gr=(x0-.5d0).*log(z1)-th.*y-x0+0.5d0.*log(2.0d0.*pi);
gi=th.*(x0-0.5d0)+y.*log(z1)-y;

for  k=1:10;
	t=z1.^(1-2.*k);
	gr=gr+a(k).*t.*cos((2.0d0.*k-1.0d0).*th);
	gi=gi-a(k).*t.*sin((2.0d0.*k-1.0d0).*th);
end;

if (x <= 7.0) ;
	gr1=0.0d0;
	gi1=0.0d0;

	for  j=0:na-1;
		gr1=gr1+.5d0.*log((x+j).^2+y.*y);
		gi1=gi1+atan(y./(x+j));
	end;

	gr=gr-gr1;
	gi=gi-gi1;
end;

if (x1 < 0.0d0) ;
	z1=sqrt(x.*x+y.*y);
	th1=atan(y./x);
	sr=-sin(pi.*x).*cosh(pi.*y);
	si=-cos(pi.*x).*sinh(pi.*y);
	z2=sqrt(sr.*sr+si.*si);
	th2=atan(si./sr);
	
	if (sr < 0.0d0) ;
		th2=pi+th2;
	end;

	gr=log(pi./(z1.*z2))-gr;
	gi=-th1-th2-gi;
	x=x1;
	y=y1;
end;

if (kf == 1) ;
	g0=exp(gr);
	gr=g0.*cos(gi);
	gi=g0.*sin(gi);
end;

return;