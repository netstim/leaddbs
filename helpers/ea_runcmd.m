function varargout = ea_runcmd(cmd, env, timeout)
% Run system command constructed using external binaries.

arguments
    cmd     {mustBeTextScalar}
    env     {mustBeText} = '' % env to be overridden, can be 'key=value' or {'key1=value1', 'key2=value2'}
    timeout {mustBeTextScalar} = '' % Execute cmd with timeout, can be 10s, 10m, 1h, etc.
end

if isempty(env)
    envOverride = '';
else
    if ischar(env)
        env = {env};
    end

    if isunix
        envOverride = ['export ' strjoin(env, ' ') ';'];
    else
        envOverride = ['set "' strjoin(env, '" & set "') '" & '];
    end
end

cmd = [envOverride, cmd]; 

if isunix
    cmd = ['bash -c "', cmd, '"'];
end

if ~isempty(timeout)
    % Add timeout
    binFolder = fileparts(mfilename('fullpath'));
    if ismac
        binPath = ea_getExec(fullfile(binFolder, 'gtimeout'), escapePath=true);
        cmd = [binPath ' ' timeout ' ' cmd];
    elseif isunix
        cmd = ['timeout ' timeout ' ' cmd];
    elseif ispc
        binPath = ea_getExec(fullfile(binFolder, 'procgov64'), escapePath=true);
        cmd = [binPath ' -t ' timeout ' ' cmd];
    end
end

if nargout == 0
    system(cmd);
elseif nargout == 1
    varargout{1} = system(cmd);
elseif nargout == 2
    [varargout{1}, varargout{2}] = system(cmd);
    varargout{2} = strip(varargout{2});
end
