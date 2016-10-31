% rmHotellingT2.m: performs repeated-measures hotelling T2 test
%
%        $Id: imGradient2D.m 2733 2013-05-13 11:47:54Z julien $
%           
%   number of variables fixed to 2, overlays must be given alternating
%   between the 2 variables. If number of variables is 2, data can be given
%   as amplitude (first) and phase (second) (in which case the
%   amplitudePhase logical variable must be set to true)

function [P,N,T2] = rmHotellingT2(varargin)

justDoCentralSlice=true;

amplitudePhase=false;

[k,kk,kkk]=size(varargin{1});
p=2;
n=round(nargin/p);

if n<=p
   mrErrorDlg('(rmHotellingT2) Sample-size (n) must be greater than the number of variables (p).');
end

if nargin~=n*p
  mrErrorDlg(['(rmHotellingT2) The number of input overlays must be a multiple of the number of variables (' num2str(p) ')']);
end

if amplitudePhase && p~=2
  mrErrorDlg('(rmHotellingT2) To pass the input as amplitude and phase, the number of variables must be 2');
end

X=nan(k,kk,kkk,n,p);
for i=1:n
  for j=1:p
    X(:,:,:,i,j) = varargin{p*(i-1)+j};
  end
end
clear('varargin');    

if amplitudePhase
  for i=1:n
    temp=X(:,:,:,i,1).*exp(1i*X(:,:,:,i,2));
    X(:,:,:,i,1)=real(temp);
    X(:,:,:,i,2)=imag(temp);
  end
end

%just do central slice
if justDoCentralSlice
  mrWarnDlg('(rmHotellingT2) Restricting analysis to central slice');
  oldKkk=kkk;
  kkk=1;
  X = X(:,:,ceil(oldKkk/2),:,:);
end

X=reshape(X,[k*kk*kkk n p]);
X = permute(X,[2 3 1]);

mu=zeros([1,p]); %expected vector is 0

m=nanmean(X); %Mean vector from data matrix X.

T2=nan(k*kk*kkk,1);
N=squeeze(sum(all(~isnan(X),2))); %number of data-points at each voxel
c=0;
hWaitBar = mrWaitBar(-inf,'(rmHotellingT2) Computing statistics');
nNonNaN = nnz(N>p);
for i=find(N>p)'
  c=c+1;
  data=X(:,:,i);
  data = data(all(~isnan(data),2),:);
  S=cov(data);  %Covariance matrix from data matrix X.
  T2(i)=N(i)*(m(:,:,i)-mu)/S*(m(:,:,i)-mu)'; %Hotelling's T-Squared statistic.
  mrWaitBar( c/nNonNaN, hWaitBar);
end
mrCloseDlg(hWaitBar);

if n >= 50 %Chi-square approximation.    
  P=1-chi2cdf(T2,p); %Probability that null Ho: is true.
else  %F approximation.
  P=nan(k*kk*kkk,1);
  F=(N-p)./((N-1)*p).*T2;  
  v1=p;  %Numerator degrees of freedom.
  v2=N-p;  %Denominator degrees of freedom.
  c=0;
  hWaitBar = mrWaitBar(-inf,'(rmHotellingT2) Computing P-value');
  for i=find(N>p)'
    c=c+1;
    P(i)=1-fcdf(F(i),v1,v2(i));  %Probability that null Ho: is true.
    mrWaitBar( c/nNonNaN, hWaitBar);
  end
  mrCloseDlg(hWaitBar);
end;

T2=reshape(T2,[k kk kkk]);
P=reshape(P,[k kk kkk]);
N=reshape(N,[k kk kkk]);

if justDoCentralSlice
    T2=repmat(T2,[1 1 oldKkk]);
    T2(:,:,[1:ceil(oldKkk/2)-1 ceil(oldKkk/2)+1:11])=NaN;
    P=repmat(P,[1 1 oldKkk]);    
    P(:,:,[1:ceil(oldKkk/2)-1 ceil(oldKkk/2)+1:11])=NaN;
    N=repmat(N,[1 1 oldKkk]);
    N(:,:,[1:ceil(oldKkk/2)-1 ceil(oldKkk/2)+1:11])=NaN;
end
