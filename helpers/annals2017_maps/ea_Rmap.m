function [R,Rperm]=ea_Rmap(varargin)
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
    [X,n]=ea_genX(varargin{1:6}); % accumulates all images into image matrix X.
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
    ea_exportmap(n,R,varargin{1:6});
    
    if exist('Rd','var') % permutation test
        sRd=sort(Rperm,1,'descend');
        delp=R<sRd(round((pthresh/2)*itercount),:);
        deln=R>sRd(round((1-(pthresh/2))*itercount),:);
        del=logical(delp.*deln);
        R(del)=nan;
        
        [pth,fn,ext]=fileparts(varargin{3});
        varargin{3}=fullfile(pth,[fn,'_sig',ext]);
        ea_exportmap(n,R,varargin{1:6});
    end

    
    