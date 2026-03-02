function wslpath = vta_windowspathstowsl(windowspath)
    % Adjust Windows paths for WSL
    % by Till Dembek
    windowspath = strrep(windowspath,'\','/');
    tmpidx = strfind(windowspath,':');
    drive = lower(windowspath(tmpidx-1));
    wslpath = ['/mnt/',drive,windowspath(tmpidx+1:end)];
end