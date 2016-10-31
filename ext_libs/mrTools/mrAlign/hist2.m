function [h,binWidths,bins] = hist2(v1,v2,nbins)
%
% function h = hist2(v1,v2,[nbins])
%
% Computes the joint histogram of two arrays, ignoring the tails of the
% two margingal distributions.
% Retuns a 2D array (nbins X nbins).
% Defaults nbins: sqrt(length(v1)/10)
% 
% djh, 3/2005

% Error if the two arrays are unequal sizes
if any(size(v1) ~= size(v2))
    error('Array sizes must be equal');
end

% Convert input arrays to column vectors
v1 = v1(:);
v2 = v2(:);

% Default nbins
if ~exist('nbins','var')
    nbins = round(sqrt(length(v1)/10));
end

% Initialize result array
h = zeros(nbins);

% Choose range, ignoring the tails of the two marginal distributions
histThresh = length(v1)/1000;
[cnt, val] = hist(v1,100);
goodVals = find(cnt>histThresh);
minVal1 = val(min(goodVals));
maxVal1 = val(max(goodVals));

[cnt, val] = hist(v2,100);
goodVals = find(cnt>histThresh);
minVal2 = val(min(goodVals));
maxVal2 = val(max(goodVals));

% Compute binwidths and bins
binWidth1 = (maxVal1 - minVal1)/nbins;
bins1 = [minVal1:binWidth1:maxVal1]';
binWidth2 = (maxVal2 - minVal2)/nbins;
bins2 = [minVal2:binWidth2:maxVal2]';
binWidths = [binWidth1 binWidth2];
bins = [bins1, bins2];

% Bin them
v1bin = round((v1 - minVal1)/binWidth1);
v2bin = round((v2 - minVal2)/binWidth2);

% Count them
for id = 1:length(v1)
    i = v1bin(id);
    j = v2bin(id);
    if ((1<=i) & (i<=nbins) & (1<=j) & (j<=nbins))
        h(i,j) = h(i,j) + 1;
    end
end

return

% Test
v1 = randn(1e5,1);
v2 = 100*randn(1e5,1) + 10;
[h,binWidths,bins] = hist2(v1,v2);
figure(1)
showIm(h);
figure(2)
plot(sum(h,1));
figure(3)
plot(sum(h,2));

