function varargout=ea_explorer_stats_2sampleweightedlinreg(varargin)
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
    varargout{1}.name="2-Sample Weighted Regression";
    varargout{1}.file = mfilename;
    varargout{1}.type = "2-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field','Sigmoid Field','VTA'};    
    return
end

% Actual test:
outcomein=repmat(outcomein',size(valsin,1),1);
group1=outcomein;

% we don't need group 2 in this test!
%group2=outcomein;
%group2(~isnan(valsin))=NaN;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);

mysyntax = 'outcome ~ 1 + condition';
ea_dispercent(1/size(valsin,1),'Calculating Weighted regression')
for i = 1:size(valsin,1)
    ea_dispercent(i/size(valsin,1))
    mytable = table;
    mytable.outcome= vertcat(group1(i,:)',group1(i,:)');
    mytable.weight = vertcat(1-valsin(i,:)',valsin(i,:)');
    mytable.condition = vertcat(zeros(size(valsin(i,:)')),ones(size(valsin(i,:)')));
    mymdl = fitlm(mytable,mysyntax,'Weights',mytable.weight);
    %% storing statitic values in maps
    psout(i) = mymdl.Coefficients.pValue(2);
    valsout(i) = mymdl.Coefficients.tStat(2);
end
ea_dispercent(i/size(valsin,1),'end')

% map outputs
varargout{1}=valsout;
varargout{2}=psout;
end