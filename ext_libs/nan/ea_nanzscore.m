function z=ea_nanzscore(varargin)

data=varargin{1};

dorobust=0;
if nargin>1
    if strcmp(varargin{2},'robust')
        dorobust=1;
    end
end

[datamean, datasd] = deal(zeros(1, size(data, 2)));
for ci = 1:size(data, 2)
    datawonan = data(:, ci); 
    datawonan = datawonan(~isnan(datawonan));
    if dorobust
        datamean(ci) = ea_robustmean(datawonan);
    else
        datamean(ci) = mean(datawonan);
    end
    datasd(ci) = std(datawonan); 
end 

z = ( data - repmat(datamean, size(data, 1), 1) ) ...
    ./ repmat(datasd, size(data, 1), 1); 