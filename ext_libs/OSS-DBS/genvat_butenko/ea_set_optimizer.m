function settings = ea_set_optimizer(options,settings)
% Set optimization parameters for OSS (Should be moved to GUI)
% By Butenko, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for electrode reconstruction and stimulation
    settings    % parameters for OSS-DBS simulation
end

% This should be set in GUI!
% in mA!
minCylindricCurrent = -3.0;
maxCylindricCurrent = 3.0;
minSegmentedCurrent = -1.5;
maxSegmentedCurrent = 1.5;
% set S to max values
[reconst, ~, ~, ~] = ea_get_reconstruction(options.subj.recon.recon);
% potentially buggy to use [1,2] for unilateral implantations
if nActiveSources(1) ~= 0 && nActiveSources(2) ~= 0
    [min_bound_per_contact, max_bound_per_contact, ~] = ea_get_currents_per_contact(minCylindricCurrent,maxCylindricCurrent, minSegmentedCurrent, maxSegmentedCurrent, reconst, [0,1], 1);
elseif nActiveSources(1) ~= 0
    [min_bound_per_contact, max_bound_per_contact, ~] = ea_get_currents_per_contact(minCylindricCurrent,maxCylindricCurrent, minSegmentedCurrent, maxSegmentedCurrent, reconst, 0, 1);
elseif nActiveSources(2) ~= 0
    [min_bound_per_contact, max_bound_per_contact, ~] = ea_get_currents_per_contact(minCylindricCurrent,maxCylindricCurrent, minSegmentedCurrent, maxSegmentedCurrent, reconst, 1, 1);
end

% 1-D current protocols are not processed correctly
current_dif = [max_bound_per_contact - min_bound_per_contact;max_bound_per_contact - min_bound_per_contact];
T = array2table(current_dif);
for contact = 1:size(current_dif,2)
    T.Properties.VariableNames(contact) = {['Contact_',num2str(contact-1)]};
end
writetable(T,[outputPaths.outputDir,filesep,'Current_protocols_0.csv'])
writetable(T,[outputPaths.outputDir,filesep,'Current_protocols_1.csv'])


% this should be defined in the GUI
settings.SSE_tract = zeros(size(options.prefs.machine.vatsettings.butenko_axonLength));
settings.CSE_tract = zeros(size(options.prefs.machine.vatsettings.butenko_axonLength));

% set some to SSE and CSE
settings.SSE_tract(4,1) = 1;
settings.CSE_tract(end,1) = 1;

% create a master dictionary that will guide
% in this specific case, we convert to A
jsonDict.netblendict.min_bound_per_contact = min_bound_per_contact .* 0.001;
jsonDict.netblendict.max_bound_per_contact = max_bound_per_contact .* 0.001;

% pass it as within a Docker
jsonDict.netblendict.ActivProfileDict = [outputPaths.outputDir,filesep,'target_profiles.json'];
jsonDict.netblendict.symptom_weights_file = [outputPaths.outputDir,filesep, 'Fixed_symptoms.json'];
jsonDict.netblendict.num_iterations = 1000;
jsonDict.netblendict.optim_alg = 'Dual Annealing';
jsonDict.netblendict.similarity_metric = 'Canberra';

% IMPORTANT: if StimSets is 0, but netblend is provided, then
% Launcher thinks it is a simple optimization
jsonDict.netblendict.Optimizer = 1;
jsonDict.netblendict.StimSets = 0;

netblend_settings_file = [outputPaths.outputDir,filesep, 'netblend_dict.json'];
jsonText = jsonencode(jsonDict);
fid = fopen(netblend_settings_file, 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

% don't use fixed weights here
jsonDict2.fixed_symptom_weights = [];
jsonText = jsonencode(jsonDict2);
jsonText(end-2:end-1) = '{}';

json_file = [outputPaths.outputDir,filesep,'Fixed_symptoms.json'];
fid = fopen(json_file, 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

end