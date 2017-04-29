function [options, presentfiles] = ea_assignpretra(options)

directory = [options.root,options.patientname,filesep];

presfiles=dir([directory,options.prefs.prenii_searchstring]);
pfcell = {presfiles.name}';

% order the anatomical images in accordance with 'prefs.prenii_order'
prenii_order = cellfun(@(x) strrep(options.prefs.prenii_searchstring,'*',x), options.prefs.prenii_order, 'UniformOutput', 0);
[~,idx] = ismember(prenii_order, pfcell);
presentfiles = pfcell([nonzeros(idx)',setdiff(1:numel(pfcell),nonzeros(idx))]);
    options.primarytemplate = 't2'; % default T2.

if isempty(presentfiles)
    warning(['No anatomy information found! Please put either ', ...
        prenii_order{1},', ',prenii_order{2},' or ',prenii_order{3},' into subject folder.']);
    return
end

% set prenii_unnormalized
options.prefs.prenii_unnormalized = presentfiles{1};

% determine primary template
if any(idx)
    options.primarytemplate = options.prefs.prenii_order{find(idx,1)};
end
