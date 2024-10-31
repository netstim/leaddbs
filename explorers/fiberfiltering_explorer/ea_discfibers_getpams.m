function [pamlist,FilesExist] = ea_discfibers_getpams(obj)
% Return list of PAM results

% always initialize with entries for mirrorred stims
% mirror flag will be checked in ea_discfibers_calcvals_pam
numPatient = length(obj.allpatients);
pamlist = cell(numPatient*2,2);
FilesExist = zeros(numPatient,2);  % no sense to check for mirrored here

disp('Construct PAM list...')

% IMPORTANT: if multiple pathways were used, fiberActivation files have been already merged
% in ea_discfibers_merge_pathways!
for sub=1:numPatient % Original VAT E-field

    [~,subj_tag,~] = fileparts(obj.allpatients{sub});
    subSimPrefix = [subj_tag, '_sim-'];

    % we load stim parameters and check for each side if there was a stimulation
    % merged fiberActivations always stored in MNI
    stimFolder = [obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(0), 'gs_', obj.M.guid];
    stimParams = [stimFolder, filesep, subj_tag, '_desc-stimparameters.mat'];

    if ~isfile(stimParams)
        FilesExist([sub,sub+numPatient], :) = ones(2);
        ea_cprintf('CmdWinWarnings', 'Stimulation parameters not found! Skip checking stimulation/vta existence.\n');
        return
    else
        S = ea_loadstimulation(stimParams);

        for side = 1:2
            % if no stimulation, we do not expect the file to exist
            % so no re-simulation is needed
            if all(S.amplitude{1,side} == 0)
                FilesExist(sub,side) = 1;
                pamlist{sub,side} = 'skip';

                % also mark mirrored
                if side == 1
                    pamlist{sub+numPatient,2} = 'skip';
                else
                    pamlist{sub+numPatient,1} = 'skip';
                end

                continue
            else
                if side == 1
                    BIDS_side = 'R';
                    % mirrored
                    pamlist{sub+numPatient,1} = fullfile(stimFolder,[obj.connectome, filesep, 'PAM', filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-L.mat']);
                else
                    BIDS_side = 'L';
                    % mirrored
                    pamlist{sub+numPatient,2} = fullfile(stimFolder,[obj.connectome, filesep, 'PAM', filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-R.mat']);
                end
                pamlist{sub,side} = fullfile(stimFolder,[obj.connectome, filesep, 'PAM', filesep, subSimPrefix, 'fiberActivation_model-ossdbs_hemi-',BIDS_side,'.mat']);
                FilesExist(sub,side)=exist(pamlist{sub,side},'file');
            end
        end
    end

end

FilesExist = double(logical(FilesExist)); % convert to zeros and ones.
