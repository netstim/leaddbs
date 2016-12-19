function [erg, errStr, oArg1, oArg2, oArg3]= spherical_harm(varargin)
%
%   function [erg, errStr, oArg1, oArg2, oArg3]= spherical_harm(commandStr, arg1, arg2, ....)
%
%       commandStr:
%           'init'  arg1= direction Vectors
%                   [arg2= maximum of order (default= 6)]
%                   [arg3= number of vertex for visualizing (default= 30)]
%                   [arg4= if set to 'double' all verteces are doubled in opposite direction]
%                   -------------------------------------------------------------------------
%                   erg= 1
%                   oArg1= All spherical harmonics until order= maxOrder
%                   oArg1= All inverse spherical harmonics until order= maxOrder
%
%           'getCoef'   arg1= Intesities
%                       arg2= order (<= maxOrder)
%                       -----------------------------------------------------------
%                       erg= coeficent of spherical harmonic
%                       oArg1= differences between spherical harmonics and Intesities (abs(sp) - inten) 
%
%           'getBestCoef'   arg1= Intesities
%                           [arg2= no of sample data (default= length(arg1)]
%                           [arg3= ['even' | 'even&odd' (default= 'even')]
%                           -----------------------------------------------------------
%                           erg= coeficent of spherical harmonic
%                           oArg1= differences between spherical harmonics and Intesities (abs(sp) - inten) 
%                           oArg2= order of selected spherical harmonics
%                           oArg3= coefiecent of all even orders until maxOrder
%
%           'getDist'       arg1= spherical harmonics intensities
%                           [arg2= vectors of directions (default= initiated directions)]
%                           -----------------------------------------------------------
%                           erg= array of distances for vectors
%
%           'visSphH'   arg1= sherical harmonic coeffizent
%                       [arg2= value for each dir. vector to code the color]
%                       [arg3= position of center]
%                       [arg4= axes handle]
%                       -----------------------------------------------------------
%                       erg= surface handle
%
%
% Bjoern W. Kreher
%  03/03
% 
%  UNIX


persistent maxOrder;
persistent inSphHarmMx inInvSphHarmMx;
persistent outSphHarmMx thetaMx phiMx;
persistent doubleOpt

defaultVertNo= 30;
defaultMaxOrder= 6;

errStr= ''; erg= [];
maxArg= 9; 

if length(varargin) < 2
    errStr= 'spherical_harm: There have to be at least two parameters';
    return;
end

if (~isempty(varargin)) && ischar(varargin{1})
    ftrStruct= varargin{1};
else
    errStr= 'spherical_harm: First param have to be a string';
    return;
end

commandStr= varargin{1};
argNo= 0;
for i= 1:maxArg
    if length(varargin) < (i + 1)
        param{i}= [];
    else
        param{i}= varargin{i + 1};
        argNo= i;
    end
end

if ~strcmp(commandStr, 'init') && isempty(maxOrder)
    errStr= 'spherical_harm: You have to init the function first';
    return;
end

if strcmp(commandStr, 'init')
    if isempty(param{2})
        maxOrder= defaultMaxOrder;
    else
        maxOrder= param{2};
    end

    if isempty(param{3})
        vertNo= defaultVertNo;
    else
        vertNo= param{3};
    end

    if ~isempty(param{4}) && strcmp(param{4}, 'double')
        doubleOpt= 1;
        dirVc= [param{1}; -param{1}];
    else
        doubleOpt= 0;
        dirVc= param{1};
    end
    
    inSphHarmMx= {};
    inInvSphHarmMx= {};
    for i= 0:maxOrder;
        inSphHarmMx{i+1}= private_build_sph_harm(dirVc, i);
        inInvSphHarmMx{i+1}= private_invert(inSphHarmMx{i+1});
    end

    [thetaMx, phiMx, outSphHarmMx]= private_init_view_sph_harm(vertNo, maxOrder);
    oArg1= inSphHarmMx;
    oArg2= inInvSphHarmMx;
    erg= 1;
elseif strcmp(commandStr, 'getCoef')
    if ~isempty(param{2})
        if param{2} > maxOrder
            errStr= 'spherical_harm: order is higher as max. order';
            return;
        end
        order= param{2};
    else
        order= maxOrder;
    end

    if doubleOpt
        dists= [param{1}; param{1}];
    else
        dists= param{1};
    end        
    
    [erg, errStr, oArg1]= private_getCoef(inInvSphHarmMx{order+1}, dists, inSphHarmMx{order+1});
elseif strcmp(commandStr, 'getBestCoef')
    % [coef, errStr, errVals, order, allCoef]= private_getBestCoef(invSphCell, sphCell, intensities, doubleOpt, sampleNo, even_odd)    
    [erg, errStr, oArg1, oArg2, oArg3]= private_getBestCoef(inInvSphHarmMx, inSphHarmMx, param{1}, doubleOpt, param{2}, param{3});    
elseif strcmp(commandStr, 'getDist')
    order= private_dof2order(length(param{1}));
    if ~isempty(param{2})
        sphHarmMx= private_build_sph_harm(param{2}, order);
    else
        sphHarmMx= inSphHarmMx{order + 1};
    end
        
    [erg, errStr]= private_getDists(sphHarmMx, param{1}, doubleOpt);
    
elseif strcmp(commandStr, 'visSphH')
    if ~isempty(param{2}) && (numel(param{2}) ==3)
        colorCoef= param{2};
    elseif ~isempty(param{2})
        [colorCoef, errStr]= private_getCoef(inInvSphHarmMx{5}, abs(param{2}), []);    
    else
        colorCoef= [];
    end
    erg= private_view_sph_harm(param{1}, outSphHarmMx, thetaMx, phiMx, param{4}, param{3}, colorCoef);
else
    errStr= strcat('spherical_harm: ''', commandStr, ''' is not a guilty command');
    return;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ergMx= private_build_sph_harm(vertVc, order)

vertNo= size(vertVc, 1);

if size(vertVc, 2) == 3
    [theta, phi]= cart2sph(vertVc(:,1), vertVc(:,2), vertVc(:,3));
elseif size(vertVc, 2) == 2
    theta= vertVc(:,1);
    phi= vertVc(:,2);
else
    ergMx= [];
    return
end

phiCos= cos(phi+.5*pi);
ergMx= zeros(vertNo, order^2 + 2*order + 1);

for j= 0:order
    baseAdr= j^2 + j + 1;
    LegP= reshape(legendre(j, phiCos), [j+1 vertNo]);
    
    m= 0:j;    
    normFak= sqrt((j + .5)*myfactorial(j - m)./(2*pi*myfactorial(j + m)));
    ergMx(:, baseAdr:(baseAdr + j))= exp(theta*m*1i) .* (ones(vertNo, 1)*normFak) .* LegP';
    ergMx(:, (baseAdr - 1):-1:(baseAdr - j))= (ones(vertNo, 1)*((-ones(1, j)).^m(2:end))) .* conj(ergMx(:, (baseAdr + 1):(baseAdr + j)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function invSphMx= private_invert(sphMx)
invSphMx= inv(sphMx'*sphMx)*sphMx';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dists, errStr]= private_getDists(sphMx, coef, doubleOpt)
errStr= '';
if doubleOpt
    len= size(sphMx, 1);
    dists= abs(sphMx(1:(len/2), :)*coef);
else
    dists= abs(sphMx*coef);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coef, errStr, errVals]= private_getCoef(iSphMx, intensities, sphMx)

errStr= ''; coef= [];

if numel(intensities) ~= size(iSphMx, 2)
    errStr= 'spherical_harm: Error wrong number of elements in itensitie';
    return
end

intensities= reshape(intensities, [numel(intensities), 1]);
coef= iSphMx * intensities;

if ~isempty(sphMx)
    errVals= intensities - abs(sphMx*coef);
else
    errVals= [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [thetaMx, phiMx, sphHarmMx]= private_init_view_sph_harm(vertNo, order)

[thetaMx, phiMx]= meshgrid(linspace(0, 2*pi, vertNo), linspace(-.5*pi, .5*pi, vertNo));

thetaAy= reshape(thetaMx, [numel(thetaMx) 1]);
phiAy= reshape(phiMx, [numel(phiMx) 1]);

%sphHarmMx= build_sph_harm([thetaAy phiAy], order);
sphHarmMx= private_build_sph_harm([thetaAy phiAy], order);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  erg= private_view_sph_harm(coef, sphHarmMx, thetaMx, phiMx, axesHd, pos, colorCoef)

erg= [];
coefNew= zeros(size(sphHarmMx, 2), 1);
coefNew(1:length(coef))= coef;
radMx= reshape(sphHarmMx*reshape(coefNew, [numel(coefNew) 1]), size(phiMx));

[xMx, yMx, zMx]= sph2cart(thetaMx, phiMx, abs(radMx));

if isempty(colorCoef)
%    colMx= sqrt(xMx.^2 + yMx.^2 + zMx.^2);
    colMx= zeros([size(xMx) 3]);
    len= sqrt(xMx.^2 + yMx.^2 + zMx.^2);
    minLen= min(min(len));
    maxLen= max(max(len));
    bright= (len - minLen)./(maxLen - minLen);
%    colMx(:, :, 2)= int8(1 - 2*acos(abs(xMx./len))/pi)
    colMx(:, :, 2)= bright.*abs(xMx./len);
    colMx(:, :, 1)= bright.*abs(yMx./len);
    colMx(:, :, 3)= bright.*abs(zMx./len);
elseif numel(colorCoef) == 3
    colMx= zeros([size(xMx) 3]);
    len= sqrt(xMx.^2 + yMx.^2 + zMx.^2);
    minLen= min(min(len));
    maxLen= max(max(len));
    bright= (len - minLen)./(maxLen - minLen);
    colorRGB= colorCoef/sqrt(sum(colorCoef.^2));
%    colMx(:, :, 2)= int8(1 - 2*acos(abs(xMx./len))/pi)
    colMx(:, :, 2)= bright.*abs(colorRGB(1));
    colMx(:, :, 1)= bright.*abs(colorRGB(2));
    colMx(:, :, 3)= bright.*abs(colorRGB(3));        
else
    coefColNew= zeros(size(sphHarmMx, 2), 1);
    coefColNew(1:length(colorCoef))= colorCoef;
    colMx= abs(reshape(sphHarmMx*reshape(coefColNew, [numel(coefColNew) 1]), size(phiMx)));    
end

if ~isempty(pos)
    xMx= xMx + pos(1);
    yMx= yMx + pos(2);
    zMx= zMx + pos(3);
end

if isempty(axesHd)
    figure;
    axesHd= axes;
    axis equal
    grid on;
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
end
erg= surface(xMx, yMx, zMx, abs(colMx), 'EdgeColor', 'none', 'FaceColor', 'texturemap', 'CDataMapping', 'direct', 'Parent', axesHd);
%erg=surface(xMx, yMx, zMx, colMx, 'Parent', axesHd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coef, errStr, errVals, order, allCoef]= private_getBestCoef(invSphCell, sphCell, intensities, doubleOpt, sampleNo, even_odd)

allCoef= {};
errStr= '';

if isempty(sampleNo)
    sampleNo= length(intensities);
end

maxOrder= length(invSphCell) - 1;

if strcmp(even_odd, 'even&odd')
    orderIdx= 2:1:maxOrder;
else
    orderIdx= 2:2:maxOrder;
end    

if doubleOpt
    intensities= [intensities; intensities];
    sampleNo= sampleNo*2;
end

order= 0;
ok= 1;
[coef, errStr, errVals]= private_getCoef(invSphCell{order + 1}, intensities, sphCell{order + 1});
allCoef{end + 1}= coef;
varErrVals= var(errVals);
dists= private_getDists(sphCell{order + 1}, coef, doubleOpt);                                            %test
varDist= var(dists);                                                                                  %test
dof= private_order2dof(order);

for orderNew = orderIdx
    dofNew = private_order2dof(orderNew);
    [coefNew, errStr, errValsNew] = private_getCoef(invSphCell{orderNew + 1}, intensities, sphCell{orderNew + 1});
    allCoef{end + 1} = coefNew;
    dists = private_getDists(sphCell{orderNew + 1}, coefNew, doubleOpt);                                            %test
    varDistNew = var(dists);                                                                                  %test
    varErrValsNew= var(errValsNew);
    f_test = ((sampleNo - 1 - dofNew)/(dofNew - dof))*(varDistNew - varDist)/mean(errValsNew.^2);   %%%f-test  nach Alexander, DC; MRM(2002), 48:331-340
%    f_test= ((sampleNo - 1 - dofNew)/(dofNew - dof))*(varErrVals - varErrValsNew)/varErrValsNew;   %%%f-test  analog wie bei multi-tensor
    if (f_test > 0) && (0.01 > local_myFdist(f_test, dofNew - dof, sampleNo - 1 - dofNew))
%    if (f_test > 0) & (0.05 > local_myFdist(f_test, dofNew - dof, sampleNo - 1 - dofNew))
        order = orderNew;    dof= dofNew;
        coef = coefNew;      errVals= errValsNew;
        varDist = varDistNew;         %test
        varErrVals = varErrValsNew;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dof, errStr]= private_order2dof(order)
errStr= '';
dof= order^2 + 2*order + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [order, errStr]= private_dof2order(dof)
errStr= '';
order= sqrt(abs(dof)) - 1;

if (dof <= 0) || (round(order) ~= order)
    order= [];
    errStr= 'invalid degree of freedom for spherical harmonics';
end

function pf = local_myFdist(f,v1,v2)
%local_myFdist( F, v1, v2) returns Q(F|v1,v2), the probability
%   of observing a value of F or greater in an
%   F-distribution with v1 and v2 degrees of freedom.
%
%   e.g., in Abramowitz & Stegun's Table 26.9, the value
%   F for which Q(F|v1,v2)=0.01 when v1=3, v2=20 is 4.94;
%   this can be verified by confirming that
%   local_myFdist( 4.94, 3, 20 ) returns a probability of 0.01, or 1%.
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
