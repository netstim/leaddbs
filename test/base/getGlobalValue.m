function value = getGlobalValue(globalName)
% Safely get the Value of a global variable.
%
% USAGE:
%    value = getGlobalValue(globalName)
%
% INPUTS:
%    globalName:    The name of the global variable to get the value for
%
% OUTPUT:
%    value:         The value of the requested global variable
%
% NOTE:
%    this function has been uploaded from
%    https://github.com/opencobra/cobratoolbox/blob/master/src/base/utilities/getGlobalValue.m
%    [e6b4efd]

    eval(['global ' globalName]);
    eval(['value = ' globalName ';']);
end