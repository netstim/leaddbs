function out = ea_connectome_bids_mapper(datasetDir, subjId, session, opts)
% Launch the NIfTI→BIDS mapping GUI for func/dwi, refresh mapping, then
% build BIDS-conform connectome derivatives for Lead-Connectome.
% Returns struct with derivatives/preprocessing paths.

if ~exist('session','var') || isempty(session)
    session = 'ses-preop';
end
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end

datasetDir = GetFullPath(datasetDir);
if ~isfolder(datasetDir)
    error('Dataset folder not found: %s', datasetDir);
end
if startsWith(subjId,'sub-')
    subjBare = regexprep(subjId,'^sub-','');
else
    subjBare = subjId;
end

% First validate mapping; if OK and runs found, skip GUI and build directly
report = ea_validate_bids_mapping(datasetDir, subjBare, session);
needsMapper = ~isempty(report.issues) || (isempty(report.funcRuns) && isempty(report.dwiRuns));

if needsMapper
    % Gather NIfTI candidates across ALL sessions for this subject (anat/func/dwi)
    % so the mapping GUI can show both preop and postop together
    subjRawDir = fullfile(datasetDir,'rawdata',['sub-',subjBare]);
    niiFiles = ea_regexpdir(subjRawDir, '\.nii(\.gz)?$', 1, 'f');
    if isempty(niiFiles)
        selDir = uigetdir(datasetDir, 'Choose a folder containing NIfTI files');
        if isnumeric(selDir)
            out = struct();
            return;
        end
        niiFiles = ea_regexpdir(selDir, '\.nii(\.gz)?$', 1, 'f');
    end

    % Open interactive mapping GUI (user assigns preop/postop, func/dwi, task/run)
    try
        ea_nifti_to_bids(niiFiles, datasetDir, ['sub-',subjBare]);
    catch ME
        ea_cprintf('CmdWinWarnings','Opening mapping GUI failed: %s\n', ME.message);
    end

    % Refresh mapping JSON
    try
        ea_genrawimagesjson(datasetDir, subjBare);
    catch
    end

    % Re-validate mapping for requested session
    report = ea_validate_bids_mapping(datasetDir, subjBare, session);
    if ~isempty(report.issues)
        msgbox(sprintf('BIDS mapping issues for sub-%s:\n- %s', subjBare, strjoin(report.issues,'\n- ')), 'BIDS Mapping Issues', 'warn');
    end
end

% Build derivatives (merged rs-fMRI and selected DWI) if mapping produced runs
if ~isempty(report.funcRuns) || ~isempty(report.dwiRuns)
    if ~(isfield(opts,'deferBuild') && opts.deferBuild)
        buildOpts = struct();
        % If multiple DWI runs exist, ask user to select ONE (no concatenation)
        if numel(report.dwiRuns) > 1
            try
                [sel, ok] = listdlg('PromptString','Multiple DWI runs found. Select one to use:', ...
                    'ListString', regexprep(report.dwiRuns, ['^', regexptranslate('escape', [datasetDir, filesep])], ''), ...
                    'SelectionMode','single');
                if ok && ~isempty(sel)
                    buildOpts.selectedDwi = report.dwiRuns(sel);
                else
                    % user cancelled selection; return without building
                    out = struct();
                    return;
                end
            catch
            end
        end
        ea_build_connectome_inputs(datasetDir, subjBare, session, buildOpts);
    end
end

% Return useful paths
out = struct();
out.derivativesDir = fullfile(datasetDir,'derivatives','leaddbs',['sub-',subjBare]);
out.preprocFunc = fullfile(out.derivativesDir,'preprocessing','func');
out.preprocDwi  = fullfile(out.derivativesDir,'preprocessing','dwi');


