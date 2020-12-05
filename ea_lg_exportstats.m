function ea_lg_exportstats(M, output)

if ischar(M)
    if exist(M,'dir') == 7 % M is the group folder
        load(fullfile(M, 'LEAD_groupanalysis.mat'));
    elseif exist(M,'file') == 2 % M is the path of the lead group file
        load(M);
    end
end

if ~exist('output', 'var')
    output = [M.root, 'ea_stats_export.mat'];
end

% Suppose homogeneous setting for all the patients
export.stimulation = M.stats(1).ea_stats.stimulation.label;

% Clinical variables
if isfield(M, 'clinical') && ~isempty(M.clinical)
    export.clinicalVars = M.clinical.vars';
    export.clinicalVarNames = M.clinical.labels';
    export.clinicalVarSelected = M.ui.clinicallist;
end

% Atlas information
atlases = dir([ea_space,'atlases']);
atlases = {atlases(cell2mat({atlases.isdir})).name};
atlases = atlases(cellfun(@(x) ~strcmp(x(1),'.'), atlases));
if isnumeric(M.ui.atlassetpopup) % Old format, atlassetpopup is numeric
    export.atlas = atlases{M.ui.atlassetpopup};
else % New format, atlassetpopup is atlas name
    export.atlas = M.ui.atlassetpopup;
end
export.atlasNames = regexprep(M.vilist', '\.nii(\.gz)?$', '');

% Parcellation and connectome information
if ~isempty(M.fclist)
    labeling = dir([ea_space,'labeling',filesep,'*.nii']);
    labeling = cellfun(@(x) {strrep(x, '.nii', '')}, {labeling.name});
    export.parcellation = labeling{M.ui.labelpopup};
    export.parcellationNames = M.fclist;
    if ~strcmp(M.ui.connectomename, 'Use none')
        export.connectome = M.ui.connectomename;
    end
end

% Atlas intersection and fiber counts
AtlasIntersectionVAT = cell(length(M.stats), 2);
nAtlasIntersectionVAT = cell(length(M.stats), 2);
AtlasIntersectionEfield = cell(length(M.stats), 2);
nAtlasIntersectionEfield = cell(length(M.stats), 2);
fibercounts = cell(length(M.stats), 2);
nfibercounts = cell(length(M.stats), 2);

for i=1:length(M.stats)
    stim = M.stats(i).ea_stats.stimulation;

    if isfield(stim, 'vat') && ~isempty(stim.vat)
        AtlasIntersectionVAT(i,:) = {stim.vat.AtlasIntersection};
        nAtlasIntersectionVAT(i,:) = {stim.vat.nAtlasIntersection};
    end

    if isfield(stim, 'efield') && ~isempty(stim.efield)
        AtlasIntersectionEfield(i,:) = {stim.efield.AtlasIntersection};
        nAtlasIntersectionEfield(i,:) = {stim.efield.nAtlasIntersection};
    end

    if isfield(stim, 'ft') && ~isempty(stim.ft)
        fibercounts = {stim.ft.fibercounts};
        nfibercounts = {stim.ft.nfibercounts};
    end
end

if isfield(stim, 'vat') && ~isempty(stim.vat)
    export.AtlasIntersectionVAT = AtlasIntersectionVAT;
    export.nAtlasIntersectionVAT = nAtlasIntersectionVAT;
end

if isfield(stim, 'efield') && ~isempty(stim.efield)
    export.AtlasIntersectionEfield = AtlasIntersectionEfield;
    export.nAtlasIntersectionEfield = nAtlasIntersectionEfield;
end

if isfield(stim, 'ft') && ~isempty(stim.ft)
    export.fibercounts = fibercounts;
    export.fibercounts = nfibercounts;
end

save(output, '-struct', 'export');
