function [rundir, setupdir] = ea_permtest_setup(root, name, bundle, permset, stratifyByGroup, runFolder)
% Named analysis folder ("setup") for permutation runs:
%   <root>/<name>/              one analysis setup, named by the user
%       README.md               settings of the setup + list of its runs
%       setup.mat               the settings every run in it must match
%       <runFolder>/            e.g. fibers/ (fiber maps) or prediction_k5/
%
% All runs of one setup use the same patients, scores, fiber settings,
% statistic and shuffles, so they can be reported together. A new run is
% only added to an existing setup when its settings match (otherwise the
% user is asked for a new name), and an existing run folder is never
% overwritten.
%
% root            - results root, e.g. <leadgroup dir>/fiberfiltering/permutation.
% name            - setup name: letters, digits, '-' and '_' (max 64 characters).
% bundle          - adapter bundle (ea_permtest_adapter_*.m).
% permset         - shuffle set used (ea_permtest_get_permset).
% stratifyByGroup - shuffles within cohort groups.
% runFolder       - subfolder for this run (must not exist yet).
%
% rundir   - <root>/<name>/<runFolder> (created).
% setupdir - <root>/<name>.

name = char(name);
if isempty(regexp(name, '^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$', 'once'))
    ea_error(sprintf(['Invalid analysis name ''%s''.\nUse letters, digits, ''-'' and ''_'' only ', ...
        '(at most 64 characters, starting with a letter or digit).'], local_esc(name)), 'title', 'Permutation testing');
end
setup = local_setup(bundle, permset, stratifyByGroup);
setupdir = fullfile(root, name);
setupfile = fullfile(setupdir, 'setup.mat');

if isfolder(setupdir)
    if ~isfile(setupfile)
        ea_error(sprintf(['A folder named ''%s'' already exists here and is not a permutation analysis:\n%s\n', ...
            'Please choose another name.'], local_esc(name), local_esc(setupdir)), 'title', 'Permutation testing');
    end
    S = load(setupfile, 'setup');
    differ = local_differences(S.setup, setup);
    if ~isempty(differ)
        ea_error(sprintf(['The analysis ''%s'' was created with different settings (%s).\n', ...
            'Runs in one analysis must use the same settings -- please choose a new name.'], ...
            local_esc(name), strjoin(differ, ', ')), 'title', 'Permutation testing');
    end
    if isfolder(fullfile(setupdir, runFolder))
        ea_error(sprintf(['The analysis ''%s'' already contains a ''%s'' run.\n', ...
            'Existing results are never overwritten -- please choose a new name.'], ...
            local_esc(name), runFolder), 'title', 'Permutation testing');
    end
else
    [ok, msg] = mkdir(setupdir);
    if ~ok
        ea_error(sprintf('Cannot create the analysis folder:\n%s\n%s', local_esc(setupdir), local_esc(msg)));
    end
    created = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<NASGU>
    save(setupfile, 'setup', 'created');
    ea_permtest_readme_append(setupdir, ['Permutation analysis ''', name, ''''], local_readme(bundle, permset, stratifyByGroup));
end

rundir = fullfile(setupdir, runFolder);
[ok, msg] = mkdir(rundir);
if ~ok
    ea_error(sprintf('Cannot create the run folder:\n%s\n%s', local_esc(rundir), local_esc(msg)));
end
ea_permtest_readme_append(setupdir, 'Run added', {sprintf('- `%s/` (see its README.md)', runFolder)});


function setup = local_setup(bundle, permset, stratifyByGroup)
% Everything that has to be identical for runs of one setup.
sel = bundle.patientselection;
setup = struct( ...
    'tool', bundle.tool, ...
    'mode', bundle.mode, ...
    'patients', {bundle.patientIDs(sel)}, ...
    'groups', bundle.group(sel), ...
    'outcomelabels', {bundle.outcomelabels}, ...
    'outcomes', bundle.outcomes(sel, :), ...
    'mirrorsides', logical(bundle.mirrorsides), ...
    'connectome', bundle.connectome, ...
    'metric', bundle.methodid, ...
    'stattest', bundle.stattest, ...
    'H0', char(string(bundle.H0)), ...
    'stimulationmodel', bundle.statsettings.stimulationmodel, ...
    'connthreshold', bundle.statsettings.connthreshold, ...
    'efieldthreshold', bundle.statsettings.efieldthreshold, ...
    'stratifyByGroup', logical(stratifyByGroup), ...
    'Nperm', permset.Nperm, ...
    'permset', permset.fingerprint);


function differ = local_differences(a, b)
fa = fieldnames(a);
fb = fieldnames(b);
differ = setxor(fa, fb)';
for f = intersect(fa, fb)'
    if ~isequaln(a.(f{1}), b.(f{1}))
        differ{end+1} = f{1}; %#ok<AGROW>
    end
end


function L = local_readme(bundle, permset, stratifyByGroup)
L = {};
L{end+1} = sprintf('- Tool: %s, explorer ID: `%s`, mode: %s', bundle.tool, bundle.ID, bundle.mode);
L{end+1} = sprintf('- Lead group folder: `%s`', bundle.groupdir);
L{end+1} = sprintf('- Patients: %d of %d selected, mirrored E-fields: %d', numel(bundle.patientselection), ...
    bundle.Npatients, bundle.mirrorsides);
L{end+1} = sprintf('- Scores: %s', strjoin(bundle.outcomelabels, ', '));
L{end+1} = sprintf('- Connectome: %s, fiber metric: %s, stimulation model: %s', bundle.connectome, ...
    bundle.methodid, bundle.statsettings.stimulationmodel);
L{end+1} = sprintf('- Statistic: %s, H0: %s, coverage threshold: %g%%, E-field threshold: %g', bundle.stattest, ...
    char(string(bundle.H0)), bundle.statsettings.connthreshold, bundle.statsettings.efieldthreshold);
if stratifyByGroup
    L{end+1} = sprintf('- Shuffles: within cohort groups (stratified), %d permutations, permset `%s`', ...
        permset.Nperm, permset.fingerprint);
else
    L{end+1} = sprintf('- Shuffles: across all selected patients (pooled), %d permutations, permset `%s`', ...
        permset.Nperm, permset.fingerprint);
end
L{end+1} = '- Every run in this folder uses exactly these settings (checked when a run is added; stored in setup.mat).';


function s = local_esc(s)
% ea_error passes its message through sprintf -- keep paths/messages literal.
s = strrep(strrep(s, '\', '\\'), '%', '%%');
