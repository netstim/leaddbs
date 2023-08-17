function z = ea_nanzscore(data, robust, dim)
% Computes zscores for input data
%
% data:   row vector, column vector, or matrix
%         If data is a matrix, z-scores are computed over columns by default
%         or by specified dim.
% robust: if set to 'robust' as second argument, mean is computed excluding
%         outliers

if isrow(data)
    dim = 2;
elseif iscolumn(data)
    dim = 1;
else % matrix
    if ~exist('dim', 'var')
        dim = 1;
    end
end

if exist('robust', 'var') && strcmpi(robust, 'robust')
    datamean = ea_robustmean(data, dim);
else
    datamean = mean(data, dim, 'omitnan');
end

datasd = std(data, 0, dim, 'omitnan');
z = (data-datamean) ./ datasd;
