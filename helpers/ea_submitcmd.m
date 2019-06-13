function [status, cmdout] = ea_submitcmd(cmd)

if ~ispc
    [status, cmdout] = system(['bash -c "', cmd, '"']);
else
    [status, cmdout] = system(cmd);
end
