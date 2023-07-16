function ext = ea_getBinExt
% Return file extension for external binaries.

if ispc
    ext = '.exe';
else
    ext = ['.', ea_getarch];
end