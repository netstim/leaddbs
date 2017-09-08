%% invPolyArcLength3 - Numerically calculate the parameter (t) of parameterized [0,1] -> R^3 polynomial
% *given* an arc length (measured from 0) i.e. find the corresponding t to an arc length
% (e.g. given in [mm] from the origin) or in other words numerically invert polyArcLength3
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2016 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function t = invPolyArcLength3(polyCoeff, arcLength)
t = nan(size(arcLength));
for i = 1:length(arcLength)
    if(arcLength(i) == 0)
        t(i) = 0;
    else
        if(arcLength(i)<0)
            warning('invPolyArcLength3: given arcLength is negative! Assuming approx. symmetry and using abs(arcLength). This might be wrong!');
            arcLength(i) = abs(arcLength(i));
        end
        t(i) = fminsearch( @(b)(abs(arcLength(i) - polyArcLength3(polyCoeff,0,b))), 0);
    end
end
end