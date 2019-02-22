%% polyval3 - evaluate an (independent component) polynomial in R3
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function val = polyval3(polyCoeffs, evalAt)
    x = polyval(((polyCoeffs(:,1))), evalAt); 
    y = polyval(((polyCoeffs(:,2))), evalAt);
    z = polyval(((polyCoeffs(:,3))), evalAt);
    val = [x y z];
end