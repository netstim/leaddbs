function [mac_address] = ea_getMac()
    % Run the ifconfig command and capture the output
    [status, cmdout] = system('ifconfig');
    
    % Check if the command was executed successfully
    if status == 0
        % Split the output into lines
        lines = strsplit(cmdout, '\n');
        
        % Initialize a flag to indicate whether en0 section is found
        en0_found = false;
        
        % Initialize the variable to hold the MAC address
        mac_address = '';
        
        % Loop through each line to find the en0 interface and its ether address
        for i = 1:length(lines)
            line = strtrim(lines{i});
            
            % Check if the line contains the en0 interface
            if contains(line, 'en0:')
                en0_found = true;
            end
            
            % If en0 section is found, look for the ether address
            if en0_found && contains(line, 'ether')
                % Extract the MAC address from the line
                tokens = strsplit(line);
                mac_address = tokens{2};
                break;  % Exit the loop once the MAC address is found
            end
        end
        
        % Display the MAC address
        if ~isempty(mac_address)
            fprintf('The MAC address for en0 is: %s\n', mac_address);
        else
            fprintf('MAC address for en0 not found.\n');
        end
    else
        % Display an error message if the command failed
        fprintf('Failed to run ifconfig command.\n');
    end

end
