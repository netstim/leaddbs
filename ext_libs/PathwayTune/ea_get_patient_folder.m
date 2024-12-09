function [PatientID, patientPath] = ea_get_patient_folder(patientlist, pt_i)

    if isempty(patientlist)
        [selpath]=uigetdir(' ','Please Choose a Patient Folder');
        if strcmp(string(selpath),'0')
            waitfor(msgbox('No Patient Folder was chosen!'));
            return;
        else
            patientPath = selpath;
        end
    else
        patientPath = patientlist{pt_i};  % return one patient
    end
    
    [~,foldername,~] = fileparts(patientPath);
    PatientID = foldername;

end