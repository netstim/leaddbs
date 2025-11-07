function [status, output] = vta_runwslcommand(cmd)
    % Run shell commands via WSL
    % by Till Dembek
    full_cmd = sprintf('wsl bash -l -c "%s"', cmd);
    [status, output] = system(full_cmd);
end