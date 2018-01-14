function [T,Mn,Medn,N,nii]=ea_getMvat(M,Xs)
X=Xs{1};
XR=Xs{2};
selectedregressor=M.clinical.vars{M.ui.clinicallist};
selectedregressor=selectedregressor(M.ui.listselect,:);
percthresh=0.1; % Voxels need to be covered by how many percent of VTAs to be taken into account.

if size(selectedregressor,2)==1
    bihemispheric=0;
elseif size(selectedregressor,2)==2
    bihemispheric=1;
else
    ea_error('Please select a regressor with entries for each hemisphere or each patient to perform this action.');
end

N=sum(X,1); % N-image
patnum=length(M.ui.listselect);

X=X.*repmat(selectedregressor(:,1),1,size(X,2)); % add regressor to VTAs
X(X==0)=nan;
[~,~,~,statT]=ttest(X);
T=statT.tstat; % T-image
Mn=ea_nanmean(X); % mean image
Medn=nanmedian(X); % median image


if bihemispheric % export second side
    NR=sum(XR,1); % N-image
    
    XR=XR.*repmat(selectedregressor(:,2),1,size(XR,2)); % add regressor to VTAs
    XR(XR==0)=nan;
    
    [~,~,~,statT]=ttest(XR);
    TR=statT.tstat; % T-image
    MnR=ea_nanmean(XR); % mean image
    MednR=nanmedian(XR); % median image
    
    Mn(isnan(Mn))=0; Mn=Mn+MnR;
    Medn(isnan(Medn))=0; Medn=Medn+MednR;
    T(isnan(T))=0; T=T+TR;
    N(isnan(N))=0; N=N+NR;
    
end


nii.dt=[16,0];
