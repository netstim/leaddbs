function [betanonp, tstats, Fnp, PValues, itration, SRWRfull, PValueF]=ea_nonparamreg(y,x)


%________________ Nonparametric Regression Estimation _____________________

% Inputs:
%   x as vectors of independent variables. (without vector of ones) and
%   y as a vector dependent variable.

% Outputs:
%   betanonp is the estimated coefficients; iteration is the number of
%   t-stats are t-students of the coefficients and Fnp stands for
%   nonparametric F-statstic.

% References:
% 1- David Birkes Yadolah Dodge, (1993). Alternative Methods of Regression,
% John Wiley & Sons, Inc.DOI:10.1002/9781118150238.

% This code is much faster than the old version(about 566 times in aregression with 5 
% explabnatory variables and 500 observations): 
% http://fmwww.bc.edu/repec/bocode/n/nonparamreg.m
% Copyright(c) Shapour Mohammmadi,University of Tehran, 2020.
% shmohmad@ut.ac.ir
%__________________________________________________________________________
ry=length(y);
[rx, cx]=size(x);
%finding first estimates by OLS method
b0=regress(y,[ones(ry,1) x]);
%defining beta star
bstar(:,1)=b0(2:cx+1,1);
%subtracting Means from independent variables 
 mu=mean(x);
for i=1:cx
Xc(:,i)=x(:,i)-mu(1,i);
end
%_____________itrations for finding Nonparamerteric estimates______________
tstar=1;
itration=0;
while tstar>.000001 && itration<500
    itration=itration+1;
    z=y-x*bstar;
    rankkk=sort(z);
    
    for i=1:ry
    IIII=find(z(i,1)==rankkk);
    Index(i,1)=mean(IIII);
    end

    u=Index(:,1)-0.5*(ry+1);
    d=(Xc'*Xc)^(-1)*Xc'*u;
    w=x*d;
  
  for s=1:ry-1
         k=s+1;
    weight1(s,1)=sum(abs(w(k:ry,1)-w(s,1))); 
  end  

    sumweight=sum(weight1);
    tstar1=[];
    weight=[];
 
    for s=1:ry-1
         k=s+1;
                wdiff=w(k:ry)-w(s,1);
                ww=w(k:ry,1);
                zz=z(k:ry,1);
                www=ww(wdiff~=0);
    tstar10=(zz(wdiff~=0)-z(s,1))./(www-w(s,1));
    weight0=abs(www-w(s,1))/sumweight;
        tstar1=[tstar1;tstar10];
        weight=[weight;weight0]; 
    end

    [sortedtstar, IX]=sort(tstar1);
     sortedweights=weight(IX);
     
    cumsumweights=cumsum(sortedweights);
    Index3=find(cumsumweights>0.5);
    tstar=sortedtstar(Index3(1,1));
    bstar=bstar+tstar*d;
end
%___________final values are Nonparametric estimates_______________________

%estimation of BETA0:for the estimation of BETA0 one should get median of
%Y-X*BETA^ where BETA^ is the estimated values of slopes(without BETA0)

Yhat1=y-x*bstar;
beta0=median(Yhat1);
betanonp=[beta0;bstar];


%_______________Estimation of t-stats and F________________________________

%t-stats
et=y-[ones(rx,1) x]*betanonp;
nt=length(et);
Aij=[];
for it=1:nt
     jt=it;
        Aij0=(et(it,1)+et(jt:nt,1))/2;
    Aij=[Aij;Aij0];
end
A=sort(Aij);
a=(nt*(nt+1))/4;
b=((nt*(nt+1)*(2*nt+1))/24)^0.5;
k1=round(0.5+a-1.645*b);
k2=round(0.5+a+1.645*b);
p=cx;
f=(nt/(nt-(p+1)))^0.5;
tauhat=(f*(nt^0.5)*(A(k2)-A(k1)))/(2*1.645);
Q=[ones(rx,1) x];
VC=(tauhat^2)*((Q'*Q)^(-1));
stderrors=(diag(VC).^0.5);
tstats=betanonp./diag(VC.^0.5);
PValues=2*(1-tcdf(abs(tstats),ry-cx));

%calculation of F
efull=y-[ones(rx,1) x]*betanonp;
    rankfull=sort(efull);
    for ifull=1:ry
    Indfull=find(efull(ifull,1)==rankfull);
    Indexfull(ifull,1)=mean(Indfull);
    end
    SRWRfull=sum((Indexfull(:,1)-0.5*(ry+1)).*efull);

ereduc=y;
    rankreduc=sort(ereduc);
    for ireduc=1:ry
    Indreduc=find(ereduc(ireduc,1)==rankreduc);
    Indexreduc(ireduc,1)=mean(Indreduc);
    end
    SRWRreduc=sum((Indexreduc(:,1)-0.5*(ry+1)).*ereduc);
    
    c=(ry+1)/(48^0.5);
    Fnp=(SRWRreduc-SRWRfull)/(cx*c*tauhat);
    PValueF=1-fcdf(Fnp,cx,ry-cx-1);

%_______________________Display REsults____________________________________

disp(' ')
disp('  Results of Nonparametric Regression       ' )
disp('___________________________________________ ')
disp('   Coef.     Std.Err.   t-stats   PValues')
disp('___________________________________________ ')
disp(  [ betanonp     ,      stderrors   ,        tstats,  PValues ] )
disp('___________________________________________ ')
disp('    Tauhat   Fnp        PValue.F')
disp([tauhat          ,        Fnp      ,          PValueF])
disp('___________________________________________ ')
disp(' ')
    
