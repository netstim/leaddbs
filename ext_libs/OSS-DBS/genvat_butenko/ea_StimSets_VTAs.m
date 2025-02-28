
% Script to generate stim.volumes for StimSets or an arbitrary stim_vector
% By Butenko, konstantinmgtu@gmail.com

e_field_file = '/home/forel/Documents/data/CologneStimFit/derivatives/leaddbs/sub-CIR02/stimulations/MNI152NLin2009bAsym/stimsetstest/sub-CIR02_sim-4D_efield_model-ossdbs_hemi-R.nii';
stim_vector = [];  % in mA!
prot_index = -1;  % -1 if StimSets, otherwise provide!
StimSetsFile = '/home/forel/Documents/data/CologneStimFit/derivatives/leaddbs/sub-CIR02/stimulations/MNI152NLin2009bAsym/stimsetstest/NB_rh/Current_protocols_0.csv';

% pre-load solution over contacts
if prot_index == -1
    % use StimSets
    StimSets = table2array(readtable(StimSetsFile));
    N_contacts = size(StimSets,2);
    N_prot = size(StimSets,1);
else
    N_contacts = length(currents);
    N_prot = 1;
end

% if checking one current protocol at a time, this cell should be passed from the
% optimizer to avoid re-loading!
efield_contact_solutions = cell(N_contacts,1);
for i = 1:N_contacts
    % hardcoded naming convention
    efield_contact_solution_file = [e_field_file(1:end-4),'_',num2str(i),'.nii'];
    efield_contact_solutions{i,1} = ea_load_nii(efield_contact_solution_file);
end

for stim_i = 1:N_prot
    if prot_index == -1
        stim_vector = StimSets(stim_i, :) * 0.001;  % convert to A
        ea_get_ossdbs_StimSets_VTA(stim_vector, efield_contact_solutions,stim_i)
    else
        stim_vector = stim_vector * 0.001;  % convert to A
        ea_get_ossdbs_StimSets_VTA(stim_vector, efield_contact_solutions,prot_index)
    end
end
