function ea_save_fixed_weights(stimfolder, opt_weights, fixed_weights_values, slider_name)

% save weights that should not be changed during optimization (fixed weights)

% opt_weights - adjustable weights
% fixed_weights_values - from the sliders in ClearTune
% slider_name - corresponds to the item (symptom) name

% create a json with fixed symptoms and weights      
for i = 1:length(fixed_weights_values) 
    if opt_weights{i}.Value == 0

        % for now, weight the same for both sides
        jsonDict.fixed_symptom_weights.(genvarname([slider_name{i},'_rh'])) = fixed_weights_values(i);
        jsonDict.fixed_symptom_weights.(genvarname([slider_name{i},'_lh'])) = fixed_weights_values(i);    
    end
end

% if no fixed symptoms, create an empty dictionary
if ~exist('jsonDict','var')
    jsonDict.fixed_symptom_weights = [];
    jsonText = jsonencode(jsonDict);
    jsonText(end-2:end-1) = '{}';
else
    jsonText = jsonencode(jsonDict);
end

json_file = [stimfolder,filesep,'Fixed_symptoms.json'];
fid = fopen(json_file, 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

end