function binPath = ea_findBinPath(bin)
% Find the path containing bin

% Command to locate binary is different depending on OS
if isunix
    cmd = 'which';
else
    cmd = 'where';
end

[status, binPath] = system([cmd, ' ', bin]);
if status % Not found
    binPath = '';
else % Found, return the path
    binPath = splitlines(strip(binPath));
    binPath = fileparts(binPath{1}); % Only return first entry
end
