function ea_set_optimizer(options,stim_folder)
% Set optimization parameters for OSS (Should be moved to GUI)
% By Butenko, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    stim_folder % stimulation folder
end

% create NB folders
if exist([stim_folder,filesep,'NB_rh'],"dir")
    ea_delete([stim_folder,filesep,'NB_rh'])
end
mkdir([stim_folder,filesep,'NB_rh'])

if exist([stim_folder,filesep,'NB_lh'],"dir")
    ea_delete([stim_folder,filesep,'NB_lh'])
end
mkdir([stim_folder,filesep,'NB_lh'])

% This is just to create a dummy file, actual optimization parameters are
% in netblend_dict. 
% Actually, these are used for pre-filtering
minCylindricCurrent = -5.0;
maxCylindricCurrent = 5.0;
minSegmentedCurrent = -5.0;
maxSegmentedCurrent = 5.0;
% set S to max values
[reconst, el_model_right, el_model_left, side] = ea_get_reconstruction(options.subj.recon.recon);
% potentially buggy to use [1,2] for unilateral implantations
if isempty(el_model_right)
    side = 1;
end   

[min_bound_per_contact, max_bound_per_contact, ~] = ea_get_currents_per_contact(minCylindricCurrent,maxCylindricCurrent, minSegmentedCurrent, maxSegmentedCurrent, reconst, side, 1);

% 1-D current protocols are not processed correctly
current_dif = [max_bound_per_contact - min_bound_per_contact;max_bound_per_contact - min_bound_per_contact];
T = array2table(current_dif);
for contact = 1:size(current_dif,2)
    T.Properties.VariableNames(contact) = {['Contact_',num2str(contact-1)]};
end
writetable(T,[stim_folder,filesep,'NB_rh',filesep,'Current_protocols_0.csv'])
writetable(T,[stim_folder,filesep,'NB_lh',filesep,'Current_protocols_1.csv'])

% % weight the same for now
% jsonDict2.fixed_symptom_weights = [];
% jsonText = jsonencode(jsonDict2);
% jsonText(end-2:end-1) = '{}';
% 
% json_file = [stim_folder,filesep,'Fixed_symptoms.json'];
% fid = fopen(json_file, 'w');
% fprintf(fid, '%s', jsonText)
% fclose(fid);


end