function modelLabel = ea_lcm_resolvevatmodel(vatDir, subPrefix, vtaType)
% Resolve the simulation model used for a Lead Mapper VAT seed.

stimParams = ea_regexpdir(vatDir, 'stimparameters\.mat$', false);
if ~isempty(stimParams)
    S = ea_loadstimulation(stimParams{1});
    modelLabel = ea_simModel2Label(S.model);
    return
end

% Compiled BIDS datasets can contain VAT images without the original
% stimulation-parameter MAT file. In that case, use the model entity that
% is encoded in the BIDS VAT filenames and require one complete pair.
namePrefix = [subPrefix, '_sim-', vtaType, '_model-'];
vatFiles = dir(fullfile(vatDir, [namePrefix, '*_hemi-*.nii']));
namePattern = ['^', regexptranslate('escape', namePrefix), ...
               '([^_]+)_hemi-([LR])\.nii$'];

modelLabels = {};
hemispheres = {};
for fileNo = 1:numel(vatFiles)
    entities = regexp(vatFiles(fileNo).name, namePattern, 'tokens', 'once');
    if ~isempty(entities)
        modelLabels{end+1} = entities{1}; %#ok<AGROW>
        hemispheres{end+1} = entities{2}; %#ok<AGROW>
    end
end

models = unique(modelLabels, 'stable');
completeModels = {};
for modelNo = 1:numel(models)
    modelHemispheres = hemispheres(strcmp(modelLabels, models{modelNo}));
    if isequal(sort(modelHemispheres), {'L', 'R'})
        completeModels{end+1} = models{modelNo}; %#ok<AGROW>
    end
end

if isempty(completeModels)
    error('LeadDBS:Mapper:MissingVATPair', ...
        ['No stimulation-parameter file or complete bilateral VAT pair was found in:\n%s\n' ...
         'Expected filenames beginning with: %s'], vatDir, namePrefix);
elseif numel(completeModels) > 1
    error('LeadDBS:Mapper:AmbiguousVATModel', ...
        ['Multiple bilateral VAT models were found in:\n%s\n' ...
         'Models: %s\nRestore the stimulation-parameter file to select the intended model.'], ...
        vatDir, strjoin(completeModels, ', '));
end

modelLabel = completeModels{1};
