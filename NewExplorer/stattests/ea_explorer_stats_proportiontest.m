function varargout=ea_explorer_stats_proportiontest(varargin)
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
    varargout{1}.name="Proportion Test";
    varargout{1}.file = mfilename;
    varargout{1}.type = "Binary-Outcome Tests";
    varargout{1}.outcometype = {'binary'};
    varargout{1}.compatibility = {'VTA'};    
    return
end

% Actual test:
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);
%%%% outcomein=round(outcomein./max(outcomein)); % debugging tool to binarize outcome
if ~all(outcomein==1|outcomein==0|isnan(outcomein))
    warning('off','backtrace')
    warning('Values other than 0 or 1 detected. Outcome is not binary. Please choose different test!')
    warning('on','backtrace')
    return
else
    numtrue = sum(outcomein==1);
    numfalse = sum(outcomein==0);
    
    outcomein=repmat(outcomein',size(valsin,1),1);
    valsin=valsin .* outcomein;
    sumtrue = sum(valsin==1,2);
    sumfalse = sum(valsin==0,2);

    for i=1:size(valsin,1)
        [~,psout(i), valsout(i)]  = ea_prop_test([sumtrue(i),sumfalse(i)],[numtrue,numfalse],1);
    end
end

% map outputs
varargout{1}=valsout;
varargout{2}=psout;
end