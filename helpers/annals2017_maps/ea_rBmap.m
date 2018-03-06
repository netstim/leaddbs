function ea_rBmap(varargin)
    % creates a correlation nifti file given a set of images and a
    % regressor ("R map" in Horn 2017 AoN)
   
    % ea_Rmap(fis,regressor,outputname,mask,sk) % sk can be 'k','s','sk'
    % for smoothing and normalization options. Mask only necessary if
    % choosing 'k' option.
    
    regressor=varargin{2};
    output=varargin{3};
    
    [X,n]=ea_genX(varargin{:}); % accumulates all images into image matrix X.
    
    %b=robustfit(X',regressor)
    warning('off');
    [b,stats]=cellfun(@robustfit,cellfun(@transpose,mat2cell(X,ones(1,size(X,1))),'un',0),repmat({regressor},1,size(X,1))','Uniformoutput',0);
    warning('on');
    roB=cell2mat(b);
    roB=roB(2:2:end);
    ea_exportmap(n,roB,varargin{:});

    
    