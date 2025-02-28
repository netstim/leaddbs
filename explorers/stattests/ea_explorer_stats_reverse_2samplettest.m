function varargout=ea_explorer_stats_reverse_2samplettest(varargin)
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
    varargout{1}.name="Reverse 2-Sample T-Test";
    varargout{1}.file = mfilename;
    varargout{1}.type = "Binary-Outcome Tests";
    varargout{1}.outcometype = {'binary'};
    varargout{1}.compatibility = {'Electric Field','Sigmoid Field'};    
    return
end
% Actual test:

ImpBinary=double(outcomein>0); % make sure variable is actually binary
% restore nans
ImpBinary(isnan(outcomein))=nan;
upSet=valsin(:,ImpBinary==1)';
downSet=valsin(:,ImpBinary==0)';

[~,ps,~,stats]=ttest2(full(upSet),full(downSet)); % Run two-sample t-test across connected / unconnected values
valsout=stats.tstat';
psout=ps;
   
% map outputs
varargout{1}=valsout;
varargout{2}=psout;
end