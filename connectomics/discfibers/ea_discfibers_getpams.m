function pamlist = ea_discfibers_getpams(obj)
% Return list of VATs

if obj.M.ui.detached
    pthprefix = [fileparts(obj.leadgroup),filesep];
else
    pthprefix = '';
end

% For multiple protocols, we should look for indexed files in the stim
% folder, but mark that they are from the same patient
numPatient = length(obj.allpatients);
pamlist = cell(numPatient,2);   % no mirroring

% here we will add the missing ones (you just need to know how much you have in total, iterate, set to zero for no match)

disp('Construct PAM list...')

% IMPORTANT: if multiple pathways were used, fiberActivation files have been already merged
% in ea_discfibers_merge_pathways!
for sub=1:numPatient % Original VAT E-field

    [~,subj_tag,~] = fileparts(obj.M.patient.list{sub});
    subSimPrefix = [subj_tag, '_sim-'];
    
    pamlist{sub,1} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_right.mat'];
    pamlist{sub,2} = [pthprefix, obj.allpatients{sub},filesep, 'stimulations',filesep,...
        ea_nt(0), 'gs_',obj.M.guid,filesep,subSimPrefix, 'fiberActivation_left.mat'];
end

end
