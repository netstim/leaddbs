function bundle = ea_permtest_adapter_disctract(obj, pathopts)
% Adapter for the fiber filtering explorer (ea_disctract). Collects
% everything the tool-agnostic permutation engine needs into one neutral
% struct ("bundle"), so the engine never touches the explorer object.
% obj is not modified, apart from the legacy-settings compatibility
% conversion that ea_method2methodid applies (the explorer applies it
% itself on every recalculation).
%
% pathopts (optional struct) - where the lead group folder is, used for
% the shuffle set (<groupdir>/permsets) and the results
% (<groupdir>/fiberfiltering/permutation/<analysis name>):
%   .fibfiltPath  the .fibfilt the explorer was loaded from (recommended
%                 when loading from the command line); groupdir is two
%                 folders above it
%   .groupdir     the lead group folder directly (overrides fibfiltPath)
% Without pathopts: obj.analysispath if that file exists, else the folder
% of obj.leadgroup.
%
% One explorer state can hold several fiber analyses, mirroring the
% "group" loop of ea_discfibers_calcstats.m:
%   'Single Tract Analysis'       - 1 analysis (the response variable)
%   'Split & Color By Subscore'   - 1 analysis per subscore; all subscores
%                                   share the same patients and fibers
%   'Split & Color By Group'      - 1 analysis per patient group, each on
%                                   that group's selected patients (own
%                                   fiber selection)
% 'Split & Color By PCA' is not supported yet.
%
% The fiber selection replicates step 1 ("prefiltering") of
% ea_discfibers_calcstats.m: fibers connected to too few (and, for VTAs,
% too many) of the analysis' patients are excluded. This step does not use
% the outcome, so the analysed fibers are identical for every permutation.
% ea_permtest_check_disctract.m verifies that the engine reproduces the
% explorer's own fiber values.
%
% Further requirements: one outcome column per analysis (no hemiscores),
% no covariates (stops with an error before anything is calculated), and a
% statistical test of the families 'Correlations', '2-Sample Tests' or
% 'Binary-Outcome Tests'. One-sample tests (null = shift, needs sign
% flipping rather than shuffling) and descriptive maps are not supported.
%
% bundle fields:
%   tool, mode, ID, leadgroup, groupdir, groupdirSource, resultsdir
%   patientIDs, Npatients, patientselection, group, mirrorsides, rngseed
%   outcomes        Npatients x K, one column per outcome (response
%                   variable, or one per subscore)
%   outcomelabels   1 x K cellstr
%   analyses(a)     .key        index of the analysis in the explorer's
%                               vals{key,side} (subscore number / group label)
%                   .label      e.g. subscore name or 'Group 2'
%                   .outcomeCol column of bundle.outcomes analysed
%                   .patientselection  patients of this analysis
%                   .gpatsel    columns of fibsval used (incl. mirrored)
%                   .fiberset   index into bundle.fibersets
%   fibersets(f)    .valsin{side}      analysed fibers x numel(gpatsel)
%                   .nonemptyidx{side} rows of fibsval{side} analysed
%                   .Nfibers{side}     all connected fibers on that side
%                   .fibcell{side}     coordinates of the analysed fibers
%   connectome, connid, methodid, statsettings, stattest, statfamily,
%   statfile, H0, Nsides
%   settings        snapshot of further explorer settings (documentation)

if ~exist('pathopts', 'var') || isempty(pathopts)
    pathopts = struct;
end

if ~isa(obj, 'ea_disctract')
    ea_error('ea_permtest_adapter_disctract expects a fiber filtering explorer object (ea_disctract).');
end
if obj.runwhite
    ea_error('Permutation testing is not available for connected-fiber (white) display mode.');
end
if ~ismember(obj.multitractmode, {'Single Tract Analysis', 'Split & Color By Subscore', 'Split & Color By Group'})
    ea_error(sprintf(['Permutation testing supports ''Single Tract Analysis'', ''Split & Color By Subscore'' ', ...
        'and ''Split & Color By Group'' (current mode: ''%s'').'], obj.multitractmode));
end
if ~isempty(obj.covars)
    if ~isempty(obj.covarlabels)
        covarnames = strjoin(cellfun(@char, obj.covarlabels, 'Uni', 0), ', ');
    else
        covarnames = sprintf('%d covariate(s)', numel(obj.covars));
    end
    covarnames = strrep(strrep(covarnames, '\', '\\'), '%', '%%'); % ea_error passes the message through sprintf
    ea_error(['Permutation testing is not designed for covariates yet (selected: ', covarnames, ').\n', ...
        'No permutations were calculated. Remove the covariates in the explorer to run permutation testing.'], ...
        'title', 'Permutation testing');
end

stattests = ea_explorer_statlist(obj.responsevar);
[is, idx] = ismember(obj.statsettings.stattest, stattests.name);
if ~is
    ea_error(['Function for test ', obj.statsettings.stattest, ' missing.']);
end
statfamily = char(stattests.type(idx));
if ~ismember(statfamily, {'Correlations', '2-Sample Tests', 'Binary-Outcome Tests'})
    ea_error(sprintf(['Permutation testing by shuffling outcomes is not valid for ''%s'' (%s). ', ...
        'Supported: Correlations, 2-Sample Tests, Binary-Outcome Tests.'], ...
        obj.statsettings.stattest, statfamily));
end

connid = ea_conn2connid(obj.connectome);
if ~isfield(obj.results, connid)
    ea_error(sprintf('No fiber filtering results for connectome ''%s'' -- calculate them in the explorer first.', obj.connectome));
end
methodid = ea_method2methodid(obj);

% Fiber values -- same selection as ea_discfibers_calcstats.m
switch obj.statsettings.stimulationmodel
    case 'Sigmoid Field'
        if obj.connectivity_type == 2
            fibsval = obj.results.(connid).('PAM_probA').fibsval;
        else
            fibsval_raw = obj.results.(connid).(methodid).fibsval;
            fibsval = fibsval_raw;
            for side = 1:size(fibsval_raw, 2)
                fibsval{1,side}(:,:) = ea_SigmoidFromEfield(fibsval_raw{1,side}(:,:));
            end
        end
    otherwise
        fibsval = cellfun(@full, obj.results.(connid).(methodid).fibsval, 'Uni', 0);
end

% Patient identifiers: patient folder names (machine-independent; the
% explorer re-roots all patients to the lead group's dataset folder, so
% folder names are unique within a lead group).
Npatients = numel(obj.allpatients);
patientIDs = cell(Npatients, 1);
for i = 1:Npatients
    [~, name, ext] = fileparts(regexprep(obj.allpatients{i}, '[\\/]+$', ''));
    patientIDs{i} = [name, ext];
end
if numel(unique(patientIDs)) < Npatients
    ea_error('Patient folder names are not unique in this lead group -- cannot build patient IDs.');
end

if isfield(obj.M, 'patient') && isfield(obj.M.patient, 'group') && numel(obj.M.patient.group) == Npatients
    group = obj.M.patient.group(:);
else
    group = ones(Npatients, 1);
end

patsel = obj.patientselection(:)';

% Outcomes and analyses -- same "group" loop as ea_discfibers_calcstats.m
switch obj.multitractmode
    case {'Single Tract Analysis', 'Split & Color By Group'}
        if size(obj.responsevar, 2) > 1
            ea_error('Hemiscore outcome (2 columns) is not yet supported for permutation testing.');
        end
        if numel(obj.responsevar) ~= Npatients
            ea_error('Outcome length does not match the number of patients in the lead group.');
        end
        outcomes = obj.responsevar(:);
        outcomelabels = {char(obj.responsevarlabel)};
    case 'Split & Color By Subscore'
        if ~isfield(obj.subscore, 'vars') || isempty(obj.subscore.vars)
            ea_error('''Split & Color By Subscore'' is selected but no subscores are defined.');
        end
        K = numel(obj.subscore.vars);
        outcomes = nan(Npatients, K);
        for k = 1:K
            v = obj.subscore.vars{k};
            if size(v, 2) > 1
                ea_error(sprintf('Subscore %d has 2 columns (hemiscore) -- not yet supported for permutation testing.', k));
            end
            if numel(v) ~= Npatients
                ea_error(sprintf('Subscore %d length does not match the number of patients in the lead group.', k));
            end
            outcomes(:, k) = v(:);
        end
        if isfield(obj.subscore, 'labels') && numel(obj.subscore.labels) == K
            outcomelabels = cellfun(@char, obj.subscore.labels(:)', 'Uni', 0);
        else
            outcomelabels = arrayfun(@(k) sprintf('Subscore %d', k), 1:K, 'Uni', 0);
        end
end

analyses = struct('key', {}, 'label', {}, 'outcomeCol', {}, 'patientselection', {}, 'gpatsel', {}, 'fiberset', {});
switch obj.multitractmode
    case 'Single Tract Analysis'
        analyses(1) = local_analysis(1, outcomelabels{1}, 1, patsel, 1, obj.mirrorsides, Npatients);
    case 'Split & Color By Subscore'
        for k = 1:size(outcomes, 2)
            analyses(k) = local_analysis(k, outcomelabels{k}, k, patsel, 1, obj.mirrorsides, Npatients);
        end
    case 'Split & Color By Group'
        if ~isfield(obj.M, 'patient') || ~isfield(obj.M.patient, 'group')
            ea_error('''Split & Color By Group'' is selected but the lead group has no patient groups.');
        end
        glabels = unique(obj.M.patient.group(obj.patientselection))';
        for gi = 1:numel(glabels)
            groupspt = find((obj.M.patient.group)' == glabels(gi));
            gsel = groupspt(ismember(groupspt, patsel));
            analyses(gi) = local_analysis(glabels(gi), sprintf('Group %d', glabels(gi)), 1, gsel, gi, obj.mirrorsides, Npatients);
        end
end

% Fiber sets: one per distinct patient selection (prefiltering depends on
% the patients, not on the outcome)
Nsides = numel(fibsval);
connthr = obj.statsettings.connthreshold / 100;
Nfibersets = max([analyses.fiberset]);
fibersets = struct('valsin', cell(1, Nfibersets), 'nonemptyidx', [], 'Nfibers', [], 'fibcell', []);
for f = 1:Nfibersets
    gpatsel = analyses(find([analyses.fiberset] == f, 1)).gpatsel;
    for side = 1:Nsides
        gfibsval = fibsval{side};
        switch obj.statsettings.stimulationmodel
            case 'VTA'
                Nmap = sum(gfibsval(:,gpatsel), 2);
            case {'Sigmoid Field', 'Electric Field'}
                Nmap = sum((gfibsval(:,gpatsel) > obj.statsettings.efieldthreshold), 2);
            otherwise
                ea_error(['Unsupported stimulation model: ', obj.statsettings.stimulationmodel]);
        end
        gfibsval(Nmap < (connthr * length(gpatsel)), gpatsel) = nan;
        if strcmp(obj.statsettings.stimulationmodel, 'VTA')
            gfibsval(Nmap > ((1 - connthr) * length(gpatsel)), gpatsel) = nan;
        end

        nonempty = sum(gfibsval(:,gpatsel), 2, 'omitnan') > 0;
        fibersets(f).nonemptyidx{side} = find(nonempty);
        fibersets(f).valsin{side} = full(gfibsval(nonempty, gpatsel));
        fibersets(f).Nfibers{side} = size(gfibsval, 1);
        fibersets(f).fibcell{side} = obj.results.(connid).fibcell{side}(fibersets(f).nonemptyidx{side});
    end
end

% Lead group folder (shuffle sets and results are saved below it). The
% paths stored inside the explorer can be outdated when the .fibfilt was
% moved or loaded from the command line, so prefer where the .fibfilt
% actually is: <groupdir>/fiberfiltering/<ID>.fibfilt.
if isfield(pathopts, 'groupdir') && ~isempty(pathopts.groupdir)
    groupdir = char(pathopts.groupdir);
    groupdirSource = 'given explicitly';
elseif isfield(pathopts, 'fibfiltPath') && ~isempty(pathopts.fibfiltPath)
    groupdir = fileparts(fileparts(char(pathopts.fibfiltPath)));
    groupdirSource = 'location of the given .fibfilt';
elseif ~isempty(obj.analysispath) && isfile(obj.analysispath)
    groupdir = fileparts(fileparts(obj.analysispath));
    groupdirSource = 'location of the explorer''s .fibfilt (obj.analysispath)';
else
    groupdir = fileparts(obj.leadgroup);
    groupdirSource = 'lead group path stored in the explorer (obj.leadgroup)';
end
if ~isfolder(groupdir)
    ea_error(sprintf(['The lead group folder (%s) does not exist on this machine:\n%s\n', ...
        'Pass the .fibfilt you loaded or the lead group folder, e.g. ', ...
        'ea_permtest_adapter_disctract(obj, struct(''fibfiltPath'', ''<path>'')).'], ...
        groupdirSource, strrep(strrep(groupdir, '\', '\\'), '%', '%%')), 'title', 'Permutation testing');
end

bundle = struct;
bundle.formatversion = 2;
bundle.tool = 'fiberfiltering_explorer';
bundle.mode = obj.multitractmode;
bundle.ID = obj.ID;
bundle.leadgroup = obj.leadgroup;
bundle.groupdir = groupdir;
bundle.groupdirSource = groupdirSource;
bundle.resultsdir = fullfile(groupdir, 'fiberfiltering', 'permutation'); % analyses: <resultsdir>/<name>/

bundle.patientIDs = patientIDs;
bundle.Npatients = Npatients;
bundle.patientselection = patsel;
bundle.group = group;
bundle.mirrorsides = logical(obj.mirrorsides);
bundle.rngseed = obj.rngseed;

bundle.outcomes = outcomes;
bundle.outcomelabels = outcomelabels;
bundle.analyses = analyses;
bundle.fibersets = fibersets;

bundle.connectome = obj.connectome;
bundle.connid = connid;
bundle.methodid = methodid;
bundle.statsettings = obj.statsettings;
bundle.stattest = obj.statsettings.stattest;
bundle.statfamily = statfamily;
bundle.statfile = char(stattests.file(idx));
bundle.H0 = obj.statsettings.H0;
bundle.Nsides = Nsides;

bundle.settings = struct( ...
    'connectivity_type', obj.connectivity_type, ...
    'e_field_metric', obj.e_field_metric, ...
    'calculationMethod', obj.calculationMethod, ...
    'native', obj.native, ...
    'multi_pathways', obj.multi_pathways);


function an = local_analysis(key, label, outcomeCol, patsel, fiberset, mirrorsides, Npatients)
patsel = patsel(:)';
gpatsel = patsel;
if mirrorsides
    gpatsel = [gpatsel, gpatsel + Npatients];
end
an = struct('key', key, 'label', label, 'outcomeCol', outcomeCol, ...
    'patientselection', patsel, 'gpatsel', gpatsel, 'fiberset', fiberset);
