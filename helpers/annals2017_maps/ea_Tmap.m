function ea_Tmap(varargin)
    % creates a simple weighted T-map.
    
    % ea_Tmap(fis,regressor,outputname,mask,sk) % sk can be 'k','s','sk'
    % for smoothing and normalization options. Mask only necessary if
    % choosing 'k' option.
    
    regressor=varargin{2};
    output=varargin{3};
    [X,n]=ea_genX(varargin{:}); % accumulates all images into image matrix X.
    
    if isempty(regressor) % create simple average map without any weighting
        regressor=ones(size(X,2),1);
    end
    regressor=regressor-min(regressor); % shift up in case of negatives
    regressor=regressor+eps; % shift up marginally so there is no zero
    nregressor=regressor./ea_nansum(regressor); % divide by sum so sum is 1
    [h,p,ci,stats]=ttest([X.*repmat(nregressor',size(X,1),1)]');
    
    ea_exportmap(n,stats.tstat,varargin{:});

    


