function ea_save_fixed_weights(stimfolder, opt_weights, fixed_weights_values, slider_name)

% save weights that should not be changed during optimization (fixed weights)
% opt_weights - adjustable weights

% create a json with fixed symptoms and weights      
for i = 1:length(app.tractset.subscore.weights) 
    if opt_weights{i}.Value == 1

        % for now, weight the same for both sides
        jsonDict.fixed_symptom_weights.(genvarname([slider_name{i},'_rh'])) = fixed_weights_values(i);
        jsonDict.fixed_symptom_weights.(genvarname([slider_name{i},'_lh'])) = fixed_weights_values(i);    
    end
end

json_file = [stimfolder,filesep,'Fixed_symptoms.json'];

jsonText = jsonencode(jsonDict);
fid = fopen(json_file, 'w');
fprintf(fid, '%s', jsonText)
fclose(fid);

end