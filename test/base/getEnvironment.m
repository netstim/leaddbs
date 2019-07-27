function environment = getEnvironment()
% Get all values of current globals in a struct.
% USAGE:
%    environment = getEnvironment()
%
% OUTPUT:
%    environment:      a struct with two fields
%                       * .globals - contains all global values
%                       * .path - contains the current path
% NOTE:
%    this function has been uploaded from
%    https://github.com/opencobra/cobratoolbox/blob/master/src/base/utilities/getEnvironment.m
%    [cadfa58]

environment = struct();
globals = struct();
globalvars = who('global');
for i = 1:numel(globalvars)
    globals.(globalvars{i}) = getGlobalValue(globalvars{i});
end
environment.globals = globals;
environment.path = path;
end

