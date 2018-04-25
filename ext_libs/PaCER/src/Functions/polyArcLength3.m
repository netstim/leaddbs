%% Calculate the arc length of an interval on a parameterized polynomial in R3
% Params:
%     polyCoeff - coefficient matrix
%     lowerLimit
%     upperLimit
%
% lowerLimit / upperLimit parameters can take scalar as well as vector inputs,
% in the vector case arc lengths are returend for all pairs of intervals
% between lowerLimit and upperLimit
%
% Andreas Husch
% Centre Hospitalier de Luxembourg / Luxembourg Centre for Systems
% Biomedicine, University of Luxembourg
% 2016  - 2017
% mail@andreashusch.de, husch.andreas@chl.lu


function arcLength = polyArcLength3(polyCoeff, lowerLimit, upperLimit)
epsilon = 0.001; % avoid numerical accuracy problems in assertion
assert(all(lowerLimit(:) <= upperLimit(:) + epsilon));

regX = polyCoeff(:,1);
regY = polyCoeff(:,2);
regZ = polyCoeff(:,3);

x_d = polyder(regX);
y_d = polyder(regY);
z_d = polyder(regZ);

arcLength = nan(size(lowerLimit));

for i = 1:length(lowerLimit)
    f = @(t) sqrt(polyval(x_d, t).^2 + polyval(y_d, t).^2 + polyval(z_d, t).^2); % The arc length is defined as the integral of the norm of the derivatives of the parameterized equations.
    arcLength(i) = integral(f,lowerLimit(i),upperLimit(i));
end

end