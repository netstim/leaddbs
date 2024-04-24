function varargout=ea_explorer_stats_1sampleweightedlinreg(varargin)
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
    varargout{1}.name="1-Sample Weighted Regression";
    varargout{1}.file = mfilename;
    varargout{1}.type = "1-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field','Sigmoid Field','VTA'};    
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

%valsout = ~valsin;
%valsin_new = ones(size(valsin));
%valsout_new = ones(size(valsin));

% zero entries are meaningless!
%valsin_new(valsin==0) = nan;
%valsout_new(valsout==0) = nan;

group1=outcomein;
group1(isnan(valsin))=NaN;
group2=repmat(H0,size(valsin));
group2(isnan(group1))=NaN;
valsout=nan(size(valsin,1),1);
psout=nan(size(valsin,1),1);

mysyntax = 'outcome ~ 1 + condition';
ea_dispercent(1/size(valsin,1),'Calculating Weighted regression')
for i = 1:size(valsin,1)
    ea_dispercent(i/size(valsin,1))
    mytable = table;
    mytable.outcome= vertcat(group2(i,:)',group1(i,:)');
    mytable.weight = vertcat(valsin(i,:)',valsin(i,:)');
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
