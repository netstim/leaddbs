function ea_run_cleartune_from_leadDBS(~,~,handles, file)
    
    config = [];

    % Choose if/which fibfilt file to run with
    if file == 1
        tractset_path=ea_regexpdir(ea_get_cleartune_root, 'PD_Symptom_Specific_Tracts_Rajamani_2024.fibfilt', 1, 'file');
        if isempty(tractset_path)
            disp('Tractset file does not exist. Select another fibfilt file within the tool to proceed.');
        else
            disp('Loading the tractset file, this may take a while.');
            config = load(tractset_path{1}, '-mat');
        end        
    end
    
    % if no patient selected, just open the default cleartune window 
    uipatdir=getappdata(handles.leadfigure,'uipatdir');

    % check if uipatdir is empty
    if isempty(uipatdir)
        disp(['No patient selected. ', ...
              'Select a patient in the cleartune tool.']);
    else
        if size(uipatdir, 2) > 1
            disp(['Only one patient can be selected for cleartune. ', ...
                  'The first from the selection will be used.']);
        end
        config.patientlist = uipatdir(1); % assign the first patient
    end

    % Get the 1st patient from the list to run cleartune on    
    % Check if the reconstruction folder exists
    if ~isempty(uipatdir)
        recon_file = ea_regexpdir(uipatdir{1}, '^sub-.*_desc-reconstruction.mat$', 1, 'file'); % checks for the reconstruction file
        if isempty(recon_file)
            errordlg('No "_desc-reconstruction.mat" file found in the reconstruction folder. Please run reconstruction first.', 'File Not Found');
            return;
        end
    end

    % Generate new stimulation folder name for the patient
    config.stimulationID = char(datetime('now'), 'yyyyMMddHHmmss');

    % run cleartune
    app = ea_cleartune(config);

end