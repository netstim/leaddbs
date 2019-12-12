function result=ea_mahalanobis_cov(x,val)

%MAHALANOBIS computes the (squared) distance of each observation in x
%   from the location estimate (locvct) of the data, 
%   relative to the shape of the data.  
%
% Required input arguments:
%                  x : data matrix (n observations in rows, p variables in columns)
%             locvct : location estimate of the data (p-dimensional vector)
%      cov or invcov : scatter estimate of the data or the inverse of the scatter estimate (pxp matrix)
%
% I/O: result=mahalanobis(x,locvct,'cov',covmat)
%   The user should only give the input arguments that have to change their default value.
%   The name of the input arguments needs to be followed by their value.
%   The order of the input arguments is of no importance.
%
% Examples:
%   result=mahalanobis(x,loc,'cov',covx)
%   result=mahalanobis(x,loc,'invcov',invcovx)
%
% Output:
%   A row vector containing the squared distances of all the observations to locvct.
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Katrien Van Driessen
% Revisions by Sabine Verboven
% Last update on 18/09/2003
%

%Initialisation
result=sum(x*pinv(diag(val)).*x,2)'; 
