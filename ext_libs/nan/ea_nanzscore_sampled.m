function z=ea_nanzscore_sampled(data,samples)
% zscore using a fraction only.
datawonan = data(~isnan(data));
if samples<1 % percent given
   samples=round(numel(datawonan)*samples);
end
ixs=round(linspace(1,numel(datawonan),samples));

datamean = mean(datawonan(ixs));
datasd = std(datawonan(ixs));
z = (data-datamean)/datasd;