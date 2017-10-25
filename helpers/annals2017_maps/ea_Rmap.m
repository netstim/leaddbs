function ea_Rmap(varargin)
    % creates a correlation nifti file given a set of images and a
    % regressor ("R map" in Horn 2017 AoN)
   
    % ea_Rmap(fis,regressor,outputname,mask,sk) % sk can be 'k','s','sk'
    % for smoothing and normalization options. Mask only necessary if
    % choosing 'k' option.
    
    regressor=varargin{2};
    output=varargin{3};
    
    [X,n]=ea_genX(varargin{:}); % accumulates all images into image matrix X.
    
    R=corr(regressor,X','type','Spearman','rows','pairwise');
    
    ea_exportmap(n,R,varargin{:});

    
    