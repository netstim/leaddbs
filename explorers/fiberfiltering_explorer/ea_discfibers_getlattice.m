function [vatlist,FilesExist] = ea_discfibers_getlattice(obj)
% Return list of VATs

numPatient = length(obj.allpatients);
vatlist = cell(numPatient,2);
FilesExist = zeros(numPatient,2);  % no sense to check for mirrored here

disp('Check for 4-D niftis...')
for sub=1:numPatient
    [~,subj_tag,~] = fileparts(obj.M.patient.list{sub});
    subSimPrefix = [subj_tag, '_sim-'];

    % we load stim parameters and check for each side if there was a stimulation
    stimFolder = [obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(obj.native), 'gs_', obj.M.guid];
    stimParams = [stimFolder, filesep, subj_tag, '_desc-stimparameters.mat'];

    if ~isfile(stimParams)
        FilesExist(sub, :) = [0 0];
        ea_cprintf('CmdWinWarnings', 'Stimulation parameters not found! Skip checking stimulation/vta existence.\n');
        return
    else
        S = ea_loadstimulation(stimParams);
        for side = 1:2
            % if no stimulation, we do not expect the file to exist
            % so no re-simulation is needed
            if all(S.amplitude{1,side} == 0)
                FilesExist(sub,side) = 1;
                vatlist{sub,side} = 'skip';
                continue
            else
                if side == 1
                    BIDS_side = 'R';
                else
                    BIDS_side = 'L';
                end
                vatlist{sub,side} = fullfile(stimFolder,[subSimPrefix,'4D_efield_model-ossdbs_hemi-',BIDS_side, '.nii']);
                FilesExist(sub,side) = exist(vatlist{sub,side},'file');
            end
    
        end
    end
end

FilesExist = double(logical(FilesExist)); % convert to zeros and ones.
