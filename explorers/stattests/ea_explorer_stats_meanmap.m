function varargout=ea_explorer_stats_meanmap(varargin)
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
    varargout{1}.name="Mean-Map";
    varargout{1}.file = mfilename;
    varargout{1}.type = "Descriptive";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field','Sigmoid Field','VTA'};
    return
end

% Actual test:
outcomein=repmat(outcomein',size(valsin,1),1);
valsin=~isnan(valsin); % valsin already only includes values above the threshold;
valsout=sum((outcomein.*valsin),2,'omitnan')./sum(valsin,2,'omitnan');
psout=zeros(size(valsout));

% map outputs
varargout{1}=valsout;
varargout{2}=psout;

end
