function [status, output] = vta_runwslcommand(cmd)
    % Run shell commands via WSL
    % by Till Dembek

    % use OSS env within WSL. ToDo: pass from prefs
    oss_env_path = '/mnt/c/Users/netstim/ossv2';  
    full_cmd = sprintf('wsl bash -l -c "source %s/bin/activate && %s"', oss_env_path, cmd);
    %full_cmd = sprintf('wsl bash -l -c "%s"', cmd);
    [status, output] = system(full_cmd);
end