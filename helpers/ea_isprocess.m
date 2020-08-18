function isprocess = ea_isprocess(type, process)
% Wrapper to check if process specified by PID or NAME exists

isprocess = 0;

switch lower(type)
    case 'pid'
        % Convert to str if PID is numeric
        if isnumeric(process)
            process = num2str(process);
        end

        % Check if process specified by PID exists
        if ispc
            if ~system(['tasklist /FI "PID eq ', process, '" | find "', process, '" > NUL'])
                isprocess = 1;
            end
        else
            if ~system(['ps -p ', process, ' > /dev/null'])
                isprocess = 1;
            end
        end
    case 'name'
        % Check if process specified by NAME exists
        if ispc
            if ~system(['tasklist /FI "IMAGENAME eq ', process, '" | find "', process, '" > NUL'])
                isprocess = 1;
            end
        else
            if ~system(['pgrep -x ', process, ' > /dev/null'])
                isprocess = 1;
            end
        end
end
