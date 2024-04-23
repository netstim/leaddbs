function pamlist = ea_discfibers_getpams_mirrored(obj)
% Return list of VATs

% this is a special version where the connectome is mirrored
% and we can just flip PAM results

PAM_mirror_enabled = 1;

if PAM_mirror_enabled == 1
    ea_warndlg("PAM mirroring is used, make sure the connecome is index-symmetric!")

    numPatient = length(obj.allpatients);
    pamlist = cell(numPatient*2,2);   % no mirroring

else
    % For multiple protocols, we should look for indexed files in the stim
    % folder, but mark that they are from the same patient
    numPatient = length(obj.allpatients);
    pamlist = cell(numPatient,2);   % no mirroring
end
% here we will add the missing ones (you just need to know how much you have in total, iterate, set to zero for no match)

% here we will add the missing ones (you just need to know how much you have in total, iterate, set to zero for no match)

disp('Construct PAM list...')

% IMPORTANT: if multiple pathways were used, fiberActivation files have been already merged
% in ea_discfibers_merge_pathways!
for sub=1:numPatient % Original VAT E-field

    [~,subj_tag,~] = fileparts(obj.M.patient.list{sub});
    subSimPrefix = [subj_tag, '_sim-'];
    
    pamlist{sub,1} = [obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-R.mat'];
    pamlist{sub,2} = [obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-L.mat'];


    % Mirrored PAM, here we just initialize with a counterpart
    % actual mirroring later in calcvals
    pamlist{sub+numPatient,1} = [obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-L.mat'];
    pamlist{sub+numPatient,2} = [obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_model-ossdbs_hemi-R.mat'];

end

end