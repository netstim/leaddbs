function [R,Rperm,Rnaned,Rpermnaned]=ea_Rmap(varargin)
% creates a correlation nifti file given a set of images and a
% regressor ("R map" in Horn 2017 AoN)

% ea_Rmap(fis,regressor,outputname,mask,sk,corrtype) % sk can be 'k','s','sk'
% for smoothing and normalization options. Mask only necessary if
% choosing 'k' option.

pthresh=0.05;
itercount=1000;
regressor=varargin{2};
output=varargin{3};

if nargin<6
    corrtype='Spearman';
else
    corrtype=varargin{6};
end

[X,n]=ea_genX(varargin{1:5}); % accumulates all images into image matrix X.
X=X';
nnanix=~isnan(nansum(X,1)).*(abs(nansum(X,1)))>0;

if nargin>6
    if strcmp(varargin{7},'permute')
        Rperm=nan(1000,size(X,2));
        regressorperm=repmat(regressor,1,itercount);
        for i=1:itercount
            regressorperm(:,i)=regressorperm(randperm(numel(regressor)),i);
        end
        Rperm(:,nnanix)=corr(regressorperm,X(:,nnanix),'type',corrtype,'rows','pairwise');
    end
end

R=corr(regressor,X,'type',corrtype,'rows','pairwise');
ea_exportmap(n,R,varargin{1:5});

if exist('Rperm','var') % permutation test
    Rpermnaned=[R;Rperm]; % for now, first entry is the unpermuted one.
    sRd=sort([R;Rperm],1,'descend');
    % delete values from Rpermnaned that are not significant (uncorrected):
    delp=Rpermnaned<...
        repmat(sRd(round((pthresh/2)*itercount),:),itercount+1,1);
    deln=Rpermnaned>...
        repmat(sRd(round((1-(pthresh/2))*itercount),:),itercount+1,1);
    del=logical(delp.*deln);
    Rpermnaned(del)=nan;

    Rnaned=Rpermnaned(1,:);
    Rpermnaned=Rpermnaned(2:end,:);

    [pth,fn,ext]=fileparts(varargin{3});
    varargin{3}=fullfile(pth,[fn,'_sig',ext]);
    ea_exportmap(n,Rnaned,varargin{1:5});
end
