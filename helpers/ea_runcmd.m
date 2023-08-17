function varargout = ea_runcmd(cmd)
% Run system command constructed using external binaries.

if isunix
    cmd = ['bash -c "', cmd, '"'];
end

if nargout == 0
    system(cmd);
elseif nargout == 1
    varargout{1} = system(cmd);
elseif nargout == 2
    [varargout{1}, varargout{2}] = system(cmd);
    varargout{2} = strip(varargout{2});
end
