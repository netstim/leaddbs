function fullPath = ea_getExec(basePath, options)
% Get the full path of external executable for the current OS
%
% On macOS, it will fallback to *.maci64 in case *.maca64 doesn't exist.

arguments
    basePath             {mustBeText}   % File path without extension
    options.escapePath   {mustBeNumericOrLogical} = false % Escape the full path or not
end

basePath = GetFullPath(basePath);

if ispc
    fullPath = [basePath '.exe'];
else
    arch = ea_getarch;
    fullPath = [basePath '.' arch];
    if strcmp(arch, 'maca64') && ~isfile(fullPath)
        fullPath = replace(fullPath, ".maca64", ".maci64");
    end
end

% Sanity check
if ~isfile(fullPath)
    ea_cprintf('CmdWinWarnings', 'Requested executable doesn''t exist:\n%s\n', fullPath);
    %return;
end

if options.escapePath
    fullPath = ea_path_helper(fullPath);
end
