function ea_kill(type, process)
% Wrapper to kill process by PID or NAME

switch lower(type)
    case 'pid'
        % Convert to str if PID is numeric
        if isnumeric(process)
            process = num2str(process);
        end

        % Kill process by PID if it exists
        if ispc
            if ~system(['tasklist /FI "PID eq ', process, '" | find "', process, '" > NUL'])
                system(['taskkill /F /PID ', process, '  >NUL 2>NUL']);
            end
        else
            if ~system(['ps -p ', process, ' > /dev/null'])
                system(['kill ', process]);
            end
        end
    case 'name'
        % Kill process by NAME if it exists
        if ispc
            if ~system(['tasklist /FI "IMAGENAME eq ', process, '" | find "', process, '" > NUL'])
                system(['taskkill /F /IM "', process, '" >NUL 2>NUL'])
            end
        else
            if ~system(['pgrep -x ', process, ' > /dev/null'])
                system(['killall ', process])
            end
        end
end
