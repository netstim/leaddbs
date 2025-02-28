function varargout = ea_explorer_stats_1sampleweightedlinreg(varargin)
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

if nargin > 0
    % Map inputs
    valsin = varargin{1};
    outcomein = varargin{2};
    H0 = varargin{3};
else
    varargout{1}.name = "1-Sample Weighted Regression";
    varargout{1}.file = mfilename;
    varargout{1}.type = "1-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field', 'Sigmoid Field'};    
    return;
end

% Check if weights are outside the required range [0, 1]
if any(valsin(:) < 0) || any(valsin(:) > 1)
    valsin = normalize(valsin, 'range', [0, 1]);
end

% Actual test:
if ischar(H0)
    switch H0
        case 'Average'
            H0 = mean(outcomein, 'all', 'omitnan');
        case 'Zero'
            H0 = 0;
    end
end
outcomein = repmat(outcomein', size(valsin, 1), 1);

group1 = outcomein;
group1(isnan(valsin)) = NaN;
group2 = repmat(H0, size(valsin));
group2(isnan(group1)) = NaN;

% Initialize output arrays
valsout = nan(size(valsin, 1), 1);
psout = nan(size(valsin, 1), 1);

% Check for Parallel Computing Toolbox
if license('test', 'Distrib_Computing_Toolbox')
    % Start parallel pool if not already started
    if isempty(gcp('nocreate'))
        parpool;
    end

    % Parallel processing
    parfor i = 1:size(valsin, 1)
        % Prepare data for regression
        Y = [group2(i, :)'; group1(i, :)'];
        X = [ones(size(Y)), [zeros(size(valsin(i, :)')); ones(size(valsin(i, :)'))]];
        W = [valsin(i, :)'; valsin(i, :)'];

        % Perform weighted least squares regression
        Wsqrt = diag(sqrt(W));
        Xw = Wsqrt * X;
        Yw = Wsqrt * Y;
        b = Xw \ Yw;
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / (size(X, 1) - size(X, 2));
        C = sigma2 * inv(Xw' * Xw);
        se = sqrt(diag(C));

        % Store results
        tStat = b(2) / se(2);  % t-statistic
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), size(X, 1) - size(X, 2))); % p-value
    end
else
    % Standard processing
    for i = 1:size(valsin, 1)
        % Prepare data for regression
        Y = [group2(i, :)'; group1(i, :)'];
        X = [ones(size(Y)), [zeros(size(valsin(i, :)')); ones(size(valsin(i, :)'))]];
        W = [valsin(i, :)'; valsin(i, :)'];

        % Perform weighted least squares regression
        Wsqrt = diag(sqrt(W));
        Xw = Wsqrt * X;
        Yw = Wsqrt * Y;
        b = Xw \ Yw;
        residuals = Yw - Xw * b;
        sigma2 = (residuals' * residuals) / (size(X, 1) - size(X, 2));
        C = sigma2 * inv(Xw' * Xw);
        se = sqrt(diag(C));

        % Store results
        tStat = b(2) / se(2);  % t-statistic
        valsout(i) = tStat;
        psout(i) = 2 * (1 - tcdf(abs(tStat), size(X, 1) - size(X, 2))); % p-value
    end
end

% Map outputs
varargout{1} = valsout;
varargout{2} = psout;
end