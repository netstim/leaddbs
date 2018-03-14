function ea_Cmap(varargin)
    % creates a map both connected and correlated ("combined" map in Horn
    % 2017 AoN).
   
    % ea_Cmap(fis,regressor,outputname,mask,sk,corrtype) % sk can be 'k','s','sk'
    % for smoothing and normalization options. Mask only necessary if
    % choosing 'k' option.
    
    regressor=varargin{2};
    output=varargin{3};
    [X,n]=ea_genX(varargin{:}); % accumulates all images into image matrix X.
        
    if nargin<6
        corrtype='Spearman';
    else
        corrtype=varargin{6};
    end
    R=corr(regressor,X','type',corrtype,'rows','pairwise');
    regressor=regressor-min(regressor); % shift up in case of negatives
    regressor=regressor+eps; % shift up marginally so there is no zero
    nregressor=regressor./ea_nansum(regressor); % divide by sum so sum is 1
    A=ea_nansum(X.*repmat(nregressor',size(X,1),1),2);
    
    bidir=(R'.*A)>0; % check for entries that are either positive or negative in both maps.
    
    A(~bidir)=nan;
    
    ea_exportmap(n,A,varargin{:});

    