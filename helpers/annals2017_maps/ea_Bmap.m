function ea_Bmap(varargin)
    % creates a regression nifti file + intercept term given a set of images and a
    % regressor
   
    % ea_Rmap(fis,regressor,outputname,mask,sk,corrtype) % sk can be 'k','s','sk'
    % for smoothing and normalization options. Mask only necessary if
    % choosing 'k' option.
    
    regressor=varargin{2};
    output=varargin{3};
    
    if nargin<6
        corrtype='Spearman';
    else
        corrtype=varargin{6};
    end
    [X,n]=ea_genX(varargin{:}); % accumulates all images into image matrix X.
    nz=sum(logical(X),2)>0.2*size(X,2); % at least 20% of images covered by values.
    nnz=sum(nz);
    warning('off');
    [b,dev]=cellfun(@glmfit,cellfun(@transpose,mat2cell(X(nz,:),ones(1,nnz)),'un',0),repmat({regressor},1,nnz)','Uniformoutput',0);
    warning('on');
    
    
    bmat=cell2mat(b);
    intercept=bmat(1:2:end); % b0
    b1=bmat(2:2:end); % b1
    
    % define error
    devmat=cell2mat(dev);
    
    EXP=nan(size(X,1),1); EXP(nz)=b1;
    ea_exportmap(n,EXP,varargin{:});
    [pth,fn,ext]=fileparts(varargin{3});
    varargin{3}=fullfile(pth,[fn,'_intercept',ext]);
    EXP(nz)=intercept;
    ea_exportmap(n,EXP,varargin{:});
        
    EXP(nz)=devmat;
    varargin{3}=fullfile(pth,[fn,'_deviation',ext]);
    ea_exportmap(n,EXP,varargin{:});

    
    
    
    
    