function varargout = ea_explorer_stats_2sampleweightedlinreg(varargin)
% Function to estimate a certain mass univariate test (see below) as used
% in explorer apps.
% - Expects valsin as V x N matrix where 
%                                   V is the number of
%                                       voxels/streamlines/etc and 
%                                   N is the number of 
%                                       E-Fields/VTAs/Lesions/etc.
% - Expects outcomein as improvement values with dimension N x 1 or 2.
% - Outputs valsout (test results) and psout (p-values of test results) as V x 1 or 2.

if nargin > 0
    % Map inputs
    valsin = varargin{1};
    outcomein = varargin{2};
else
    varargout{1}.name = "2-Sample Weighted Regression";
    varargout{1}.file = mfilename;
    varargout{1}.type = "2-Sample Tests";
    varargout{1}.outcometype = {'gradual'};
    varargout{1}.compatibility = {'Electric Field', 'Sigmoid Field'};    
    return;
end

% Check if weights are outside the required range [0, 1]
if any(valsin(:) < 0) || any(valsin(:) > 1)
    valsin = normalize(valsin, 'range', [0, 1]);
end


% Actual test:
outcomein = repmat(outcomein', size(valsin, 1), 1);
group1 = outcomein;

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
        Y = [group1(i, :)'; group1(i, :)'];
        X = [ones(size(Y)), [zeros(size(valsin(i, :)')); ones(size(valsin(i, :)'))]];
        W = [1 - valsin(i, :)'; valsin(i, :)'];

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
        Y = [group1(i, :)'; group1(i, :)'];
        X = [ones(size(Y)), [zeros(size(valsin(i, :)')); ones(size(valsin(i, :)'))]];
        W = [1 - valsin(i, :)'; valsin(i, :)'];

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