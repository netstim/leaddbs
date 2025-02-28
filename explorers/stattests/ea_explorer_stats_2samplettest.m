function varargout=ea_explorer_stats_2samplettest(varargin)
% Function to estimate a certain mass univariate test (see below) as used
% in explorer apps.
% - Expects valsin as V x N matrix where 
%                                   V is the number of
%                                       voxels/streamlines/etc and 
%                                   N is the number of 
%                                       E-Fields/VTAs/Lesions/etc.
% - Expects outcomein as improvement values with dimension N x 1 or 2.
% - Outputs valsout (test results) and psout (p-values of test results) as V x 1 or 2.

if nargin>0
    % map inputs
    valsin=varargin{1};
    outcomein=varargin{2};
else
    varargout{1}.name="2-Sample T-Test";
    varargout{1}.file = mfilename;
    varargout{1}.type = "2-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'VTA'};    
    return
end

% Actual test:

outcomein=repmat(outcomein',size(valsin,1),1);

valsout = ~valsin;
valsin_new = ones(size(valsin));
valsout_new = ones(size(valsin));

% zero entries are meaningless!
valsin_new(valsin==0) = nan;
valsout_new(valsout==0) = nan;

group1=valsin_new .* outcomein;
group2=valsout_new .* outcomein;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);
[~,psout,~,stats]=ttest2(group1',group2');
psout=psout';
valsout = stats.tstat';

% map outputs
varargout{1}=valsout;
varargout{2}=psout;
end