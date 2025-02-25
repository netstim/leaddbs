function ea_querytutor(bids)
    machineId = ea_getMachineId;
    jsonData = ea_getTutorData(bids); 
    baseUrl = 'http://tutor.lead-dbs.org/?userid=';
    url = 'http://tutor.lead-dbs.org:8000/my-endpoint';
    options = weboptions('RequestMethod', 'post', 'MediaType', 'application/json', 'HeaderFields', {'User-ID', machineId});
    webwrite(url, jsonData, options);
    userUrl = strcat(baseUrl, machineId);
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
