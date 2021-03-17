function z = ea_nanzscore_sampled(data, NumSamples)
% zscore samples from input data excluding NaNs

% Data without NaNs
datawonan = data(~isnan(data));
if isempty(datawonan)
    warning('Input data are all NaN');
    % Return the original input
    z = data;
    return;
end

% Sampling
if NumSamples<1 % percent given
   NumSamples = round(numel(datawonan)*NumSamples);
end
datawonan = datawonan(round(linspace(1,numel(datawonan),NumSamples)));

% zscore
z = (data-mean(datawonan))/std(datawonan);
