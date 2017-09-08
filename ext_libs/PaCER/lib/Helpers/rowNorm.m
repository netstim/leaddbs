%% rowNorm - row-wise eucliadian norm
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function d = rowNorm(A)
d = sqrt(sum(A.*A, 2));
