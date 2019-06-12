function clearGlobal(globalName)
% Safely clear a global variable.
%
% USAGE:
%    clearGlobal(globalName)
%
% INPUTS:
%    globalName:    The name of the global variable to clear.
%
% NOTE:
%    this function has been uploaded from
%    https://github.com/opencobra/cobratoolbox/blob/master/src/base/utilities/clearGlobal.m
%    [e6b4efd]

    clearvars('-global',globalName);
end