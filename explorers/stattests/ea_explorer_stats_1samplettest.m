function varargout=ea_explorer_stats_1samplettest(varargin)
% Function to estimate a certain mass univariate test (see below) as used
% in explorer apps.
% - Expects valsin as V x N matrix where 
%                                   V is the number of
%                                       voxels/streamlines/etc and 
%                                   N is the number of 
%                                       E-Fields/VTAs/Lesions/etc.
% - Expects outcomein as improvement values with dimension N x 1 or 2.
% - Expects H0 as the value to test against {'Zero','Average' or a number}
% - Outputs valsout (test results) and psout (p-values of test results) as V x 1 or 2.

if nargin>0
    % map inputs
    valsin=varargin{1};
    outcomein=varargin{2};
    H0=varargin{3};
else
    varargout{1}.name="1-Sample T-Test";
    varargout{1}.file = mfilename;
    varargout{1}.type = "1-Sample Tests";
    varargout{1}.outcometype = {'gradual','binary'};
    varargout{1}.compatibility = {'VTA'};    
    return
end

% Actual test:

if ischar(H0)
    switch H0
        case 'Average'
            H0=mean(outcomein,'all','omitnan');
        case 'Zero'
            H0=0;
    end
end
outcomein=repmat(outcomein',size(valsin,1),1);

valsout = ~valsin;
valsin_new = ones(size(valsin));
% zero entries are meaningless!
valsin_new(valsin==0) = nan;

group1=valsin_new .* outcomein;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);
[~,psout,~,stats]=ttest(group1',H0);
psout=psout';
valsout = stats.tstat';

% map outputs
varargout{1}=valsout;
varargout{2}=psout;
end
