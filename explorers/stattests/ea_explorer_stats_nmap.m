function varargout=ea_explorer_stats_nmap(varargin)
% Function to estimate a certain mass univariate test (see below) as used
% in explorer apps.
% - Expects valsin as V x N matrix where
%                                   V is the number of
%                                       voxels/streamlines/etc and
%                                   N is the number of
%                                       E-Fields/VTAs/Lesions/etc.
% - Outputs valsout (test results) and psout (p-values of test results) as V x 1 or 2.

if nargin>0
    % map inputs
    valsin=varargin{1};
else
    varargout{1}.name="N-Map";
    varargout{1}.file = mfilename;
    varargout{1}.type = "Descriptive";
    varargout{1}.outcometype = {'gradual','binary'};
    varargout{1}.compatibility = {'Electric Field','Sigmoid Field','VTA'};
    return
end

% Actual test:
valsout=sum(valsin,2,'omitnan');
psout=nan(size(valsout));
% map outputs
varargout{1}=valsout;
varargout{2}=psout;

end
