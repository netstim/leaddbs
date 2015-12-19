function z=ea_nanzscore(varargin)

data=varargin{1};
datawonan = data(~isnan(data));
dorobust=0;
if nargin>1
    if strcmp(varargin{2},'robust');
        dorobust=1;
    end
end
if dorobust
    datamean = ea_robustmean(datawonan);
else
    datamean = mean(datawonan);
end
datasd = std(datawonan);
z = (data-datamean)/datasd;