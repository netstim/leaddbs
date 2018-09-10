function output=ea_getmodality(input)
if isstruct(input) % options supplied
    directory=[input.root,input.patientname,filesep];
    options=input;
else
    directory=input;
end
if ~strcmp(directory(end),filesep)
    directory=[directory,filesep];
end
modality = ea_checkctmrpresent(directory);
modality = find(modality);
if isempty(modality)    % no postop image present
    modality = 1;    % set to MR to work it around
elseif length(modality) == 2    % both MR and CT image present
    prefs=ea_prefs;
    modality = prefs.preferMRCT;  % set the modality according to 'prefs.preferMRCT'
end

if exist('options','var') % originally, options were supplied
    
    options.modality=modality;
    output=options;
    
else
    output=modality;
end