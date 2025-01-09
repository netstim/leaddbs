function varargout=ea_explorer_stats_ranksumtest(varargin)
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
    varargout{1}.name="Wilcoxon Rank-Sum Test";
    varargout{1}.file = mfilename;
    varargout{1}.type = "2-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'VTA'};    
    return
end
% Actual test:
outcomein=repmat(outcomein',size(valsin,1),1);

psout=nan(size(valsin,1),1);
for i=1:size(valsin,1)
    group1 = outcomein(i,valsin(i,:)==1);
    group2 = outcomein(i,valsin(i,:)==0);
    [psout(i),~,stats(i)]=ranksum(group1,group2);
end
zval=[stats(:).zval]';
valsout = sign(zval).*-log10(psout);

% map outputs
varargout{1}=valsout;
varargout{2}=psout;
end
