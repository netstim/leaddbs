function binPath = ea_findBinPath(bin)
% Find the path containing bin

% Command to locate binary is different depending on OS
if isunix
    cmd = 'which';
else
    cmd = 'where';
end

[status, binPath] = system([cmd, ' ', bin]);
if status
    binPath = ''; % Not found
else
    binPath = fileparts(strip(binPath)); % Found, return the path
end
