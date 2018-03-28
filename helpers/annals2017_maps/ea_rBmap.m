function ea_rBmap(varargin)
    % creates a correlation nifti file given a set of images and a
    % regressor ("R map" in Horn 2017 AoN)
   
    % ea_Rmap(fis,regressor,outputname,mask,sk) % sk can be 'k','s','sk'
    % for smoothing and normalization options. Mask only necessary if
    % choosing 'k' option.
    
    regressor=varargin{2};
    output=varargin{3};
    
    [X,n]=ea_genX(varargin{:}); % accumulates all images into image matrix X.
    nz=sum(logical(X),2)>0.2*size(X,2); % at least 20% of images covered by values.
    nnz=sum(nz);
    
    %b=robustfit(X',regressor)
    warning('off');
    [b]=cellfun(@robustfit,cellfun(@transpose,mat2cell(X(nz,:),ones(1,nnz)),'un',0),repmat({regressor},1,nnz)','Uniformoutput',0);
    warning('on');
    roB=cell2mat(b);
    roB=roB(2:2:end);
    intercept=roB(1:2:end);
    
    EXP=nan(size(X,1),1); EXP(nz)=roB;
    ea_exportmap(n,EXP,varargin{:});
    [pth,fn,ext]=fileparts(varargin{3});
    
    EXP(nz)=intercept;
    varargin{3}=fullfile(pth,[fn,'_intercept',ext]);
    ea_exportmap(n,EXP,varargin{:});

    
    