function [jsonData] = ea_getTutorData(handles, bids)
    options = ea_handles2options(handles);
    baseFolder = strcat(bids.datasetDir, '/derivatives/leaddbs/');
    for i = 1: size(bids.subjId, 1)
        patientName = strcat('sub-', cell2mat(bids.subjId(i)));
        recoFileName = strcat(patientName, '_desc-reconstruction.mat');
        recoFileDir = strcat(baseFolder, recoFileName);
        
        % Load data from the file
            try 
                correctReco = load(recoFileDir);
                groundMarkers = correctReco.reco.mni.markers;
                            % Store data in the struct
                groundStruct.(cell2mat(patientName)).leftElectrode.head = groundMarkers(1).head;
                groundStruct.(cell2mat(patientName)).leftElectrode.tail = groundMarkers(1).tail;
                groundStruct.(cell2mat(patientName)).leftElectrode.x = groundMarkers(1).x;
                groundStruct.(cell2mat(patientName)).rightElectrode.head = groundMarkers(2).head;
                groundStruct.(cell2mat(patientName)).rightElectrode.tail = groundMarkers(2).tail;
                groundStruct.(cell2mat(patientName)).rightElectrode.x = groundMarkers(2).x;
            catch
                errordlg('Patient ', patientName, ' localization is not complete. Please complete all localizations before submitting.');
            end
    

    end

    jsonData = jsonencode(groundStruct);

end

    