function [status, output] = vta_runwslcommand(oss_env_path,cmd)
% Run shell commands via WSL
% by Till Dembek

arguments
    oss_env_path     % (Linux-format) path to the OSS conda environment
    cmd % command to run in the terminal
end

full_cmd = sprintf('wsl bash -l -c "source %s/bin/activate && %s"', oss_env_path, cmd);
%full_cmd = sprintf('wsl bash -l -c "%s"', cmd);
[status, output] = system(full_cmd);
end