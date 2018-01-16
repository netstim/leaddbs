function [T,Mn,Medn,N]=ea_getMvat(M,Xs)
X=Xs{1};
XR=Xs{2};
clear Xs
selectedregressor=M.clinical.vars{M.ui.clinicallist};
selectedregressor=selectedregressor(M.ui.listselect,:);
X=X(M.ui.listselect,:);
if ~isnan(XR)
   XR=XR(M.ui.listselect,:); 
end

if size(selectedregressor,2)==1
    bihemispheric=0;
elseif size(selectedregressor,2)==2
    bihemispheric=1;
else
    ea_error('Please select a regressor with entries for each hemisphere or each patient to perform this action.');
end

% init vars
N=zeros(1,size(X,2));
T=zeros(1,size(X,2));
Mn=zeros(1,size(X,2));
Medn=zeros(1,size(X,2));

for lat=1:1+bihemispheric
    if lat==2
        X=XR;
    end
    Ns=sum(X,1); % N-image
    nzeros=Ns>0;
    X=X(:,nzeros);
    X=X.*repmat(selectedregressor(:,lat),1,size(X,2)); % add regressor to VTAs
    X(X==0)=nan;
    
    [~,~,~,statT]=ttest(X);
    T(nzeros)=T(nzeros)+statT.tstat; % T-image
    Mn(nzeros)=Mn(nzeros)+ea_nanmean(X); % mean image
    Medn(nzeros)=Medn(nzeros)+nanmedian(X); % median image
    N=N+Ns;
end


