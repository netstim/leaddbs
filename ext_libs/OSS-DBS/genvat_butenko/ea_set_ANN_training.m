function ea_set_ANN_training(reco_file,stim_folder,minCylindricCurrent,maxCylindricCurrent,minSegmentedCurrent,maxSegmentedCurrent)
% Set ANN training parameters and prepare OSS StimSets (Should be moved to GUI)
% By Butenko, konstantinmgtu@gmail.com

arguments
    reco_file                   % Lead-DBS options for electrode reconstruction and stimulation
    stim_folder                 % stimulation folder
    minCylindricCurrent         % lower boundary for current through cylindric contacts in mA (this is also max cathodic current!)
    maxCylindricCurrent         % upped boundary (also max anodic!)
    minSegmentedCurrent
    maxSegmentedCurrent
end

%% First we generate training and test datasets (as StimSet) based on the electrode model and provided current limits

% set S to max values
[~, el_model_right, el_model_left, side] = ea_get_reconstruction(reco_file);

% we will choose number of samples depending on the electrode model
if isempty(el_model_right)
    side = 1;  % left electrode only
end

% call python script to generate Current_protocols
if length(side) == 2
    disp("Electrodes in both hemispheres")
    for side = 0:1
        system(['python', ' ', ea_getearoot, 'cleartune/PathwayTune/TrainTest_Generator.py ', ...
            stim_folder, ' ', el_model_right, ' ', num2str(side), ' ', num2str(minCylindricCurrent), ' ', num2str(maxCylindricCurrent), ' ', num2str(minSegmentedCurrent), ' ', num2str(maxSegmentedCurrent)]);	% 0 is right side, 1 is the left side here
    end
elseif side == 0
    disp("Only right electrode")
    system(['python', ' ', ea_getearoot, 'cleartune/PathwayTune/TrainTest_Generator.py ', ...
        stim_folder, ' ', el_model_right, ' ', num2str(side), ' ', num2str(minCylindricCurrent), ' ', num2str(maxCylindricCurrent), ' ', num2str(minSegmentedCurrent), ' ', num2str(maxSegmentedCurrent)]);	% 0 is right side, 1 is the left side here
elseif side == 1
    disp("Only left electrode")
    system(['python', ' ', ea_getearoot, 'cleartune/PathwayTune/TrainTest_Generator.py ', ...
        stim_folder, ' ', el_model_left, ' ', num2str(side), ' ', num2str(minCylindricCurrent), ' ', num2str(maxCylindricCurrent), ' ', num2str(minSegmentedCurrent), ' ', num2str(maxSegmentedCurrent)]);	% 0 is right side, 1 is the left side here
else
    disp("No electrodes found, check the reconstruction file")
end

% weight the same for now
jsonDict2.fixed_symptom_weights = [];
jsonText = jsonencode(jsonDict2);
jsonText(end-2:end-1) = '{}';

json_file = [stim_folder,filesep,'Fixed_symptoms.json'];
fid = fopen(json_file, 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

% OSS will run StimSets for all these protocols

end