function jsonData = ea_getTutorData(handles, bids)
    options = ea_handles2options(handles);
    baseFolder = strcat(bids.datasetDir, '/derivatives/leaddbs/');
    patientIDs = {'pt_15454', 'pt_29781', 'pt_33544', 'pt_39468', 'pt_57245', 'pt_76325', 'pt_78754', 'pt_80206', 'pt_84257', 'pt_93127'};
    for i = 1: size(bids.subjId, 1)
        patientName = strcat('sub-', cell2mat(bids.subjId(i)));
        recoFileName = strcat(patientName, '_desc-reconstruction.mat');
        recoFileDir = strcat(baseFolder, patientName, '/reconstruction/', recoFileName);
        
        % Load data from the file
            try 
                correctReco = load(recoFileDir);
                groundMarkers = correctReco.reco.mni.markers;
                            % Store data in the struct
                groundStruct.(cell2mat(patientIDs(i))).leftElectrode.head = groundMarkers(1).head;
                groundStruct.(cell2mat(patientIDs(i))).leftElectrode.tail = groundMarkers(1).tail;
                groundStruct.(cell2mat(patientIDs(i))).leftElectrode.x = groundMarkers(1).x;
                groundStruct.(cell2mat(patientIDs(i))).rightElectrode.head = groundMarkers(2).head;
                groundStruct.(cell2mat(patientIDs(i))).rightElectrode.tail = groundMarkers(2).tail;
                groundStruct.(cell2mat(patientIDs(i))).rightElectrode.x = groundMarkers(2).x;
            catch
                errordlg('Patient ', cell2mat(patientIDs(i)), ' localization is not complete. Please complete all localizations before submitting.');
            end
    

    end

    jsonData = jsonencode(groundStruct);

end

    