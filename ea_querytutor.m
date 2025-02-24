function ea_querytutor(bids)
    mac_address = ea_getMac;
    jsonData = ea_getTutorData(bids); 
    baseUrl = 'http://tutor.lead-dbs.org/?userid=';
    url = 'http://tutor.lead-dbs.org:8000/my-endpoint';
    options = weboptions('RequestMethod', 'post', 'MediaType', 'application/json', 'HeaderFields', {'User-ID', mac_address});
    webwrite(url, jsonData, options);
    userUrl = strcat(baseUrl, mac_address);
    web(userUrl, '-browser');
end


function jsonData = ea_getTutorData(bids)  
    baseFolder = fullfile(bids.datasetDir, 'derivatives', 'leaddbs');
    patientIDs = {'pt_15454', 'pt_29781', 'pt_33544', 'pt_39468', 'pt_57245', ...
                  'pt_76325', 'pt_78754', 'pt_80206', 'pt_84257', 'pt_93127'};

    groundStruct = struct();

    for i = 1:numel(bids.subjId)
        patientName = ['sub-', bids.subjId{i}];
        recoFileName = [patientName, '_desc-reconstruction.mat'];  
        recoFileDir = fullfile(baseFolder, patientName, 'reconstruction', recoFileName);
        
        try
            correctReco = load(recoFileDir);
            groundMarkers = correctReco.reco.mni.markers;

            groundStruct.(patientIDs{i}).leftElectrode = struct( ...
                'head', groundMarkers(1).head, ...
                'tail', groundMarkers(1).tail, ...
                'x', groundMarkers(1).x);
            
            groundStruct.(patientIDs{i}).rightElectrode = struct( ...
                'head', groundMarkers(2).head, ...
                'tail', groundMarkers(2).tail, ...
                'x', groundMarkers(2).x);
            
        catch
            ea_cprintf('CmdWinWarnings', 'Localization for patient %s is not complete. Skipping.\n', patientIDs{i});
            continue;
        end
    end

    jsonData = jsonencode(groundStruct);
end


function mac_address = ea_getMac
    if ispc  % Windows
        [status, cmdout] = system('getmac');
    else  % macOS & Linux
        [status, cmdout] = system('ifconfig');
    end

    mac_address = '';

    if status == 0
        if ispc  % Parse Windows output
            lines = strsplit(cmdout, '\n');
            for i = 1:length(lines)
                line = strtrim(lines{i});
                if contains(line, '-')
                    mac_address = extractBefore(line, ' '); % First part before space is MAC
                    break;
                end
            end
        else  % Parse macOS/Linux output
            lines = strsplit(cmdout, '\n');
            found_interface = false;
            
            for i = 1:length(lines)
                line = strtrim(lines{i});
                if contains(line, 'en0:') || contains(line, 'eth0:')  % Common names: en0 (Mac) & eth0 (Linux)
                    found_interface = true;
                end
                if found_interface && contains(line, 'ether')
                    tokens = strsplit(line);
                    mac_address = tokens{2};  % Extract MAC address
                    break;
                end
            end
        end

        % Output MAC address
        if ~isempty(mac_address)
            fprintf('The MAC address is: %s\n', mac_address);
        else
            fprintf('MAC address not found.\n');
        end
    else
        fprintf('Failed to retrieve MAC address.\n');
    end
end
