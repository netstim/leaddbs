function out = ea_permtest_gui_disctract(action, app, answers)
% "Permutation testing" menu of the fiber filtering explorer
% (ea_discfiberexplorer). The app only calls
%   ea_permtest_gui_disctract('install', app)
% from its startup function; everything else lives here and calls the
% tested command-line functions of explorers/permutation_testing on the
% explorer's current state (app.tractset).
%
% Menu items (action):
%   'check' Check settings / preview  - adapter + consistency check,
%                                       fibers per analysis/side, expected
%                                       disk size of a run
%   'run'   Run permutation test      - dialog (analysis name, permutations,
%                                       stratify, workers) -> ea_permtest_fibers_disctract
%   'stats' Compute statistics        - dialog (run, alpha, direction,
%                                       step-down, GPD) -> ea_permtest_fibers_stats
%   'predict' Prediction permutation test (k-fold CV) - dialog (analysis
%                                       name, permutations, stratify,
%                                       combined target) -> ea_permtest_cv_disctract;
%                                       k, iterations and post-hoc correction
%                                       are taken from the explorer's CV panel
%   'show'  Show significant fibers   - draw an exported fiber set in the
%                                       explorer's 3D figure, coloured by R
%   'hide'  Hide significant fibers
%
% answers - optional struct that replaces the dialogs (scripting/tests):
%   run:   .name, .Nperm, .stratifyByGroup, .maxWorkers
%   stats: .rundir, .alpha, .direction, .stepdown, .useTailApprox
%   predict: .name, .Nperm, .stratifyByGroup, .targetLabel ('' = none),
%            .kIter, .posthoccorrect
%   show:  .file

if ~exist('answers', 'var')
    answers = [];
end
out = [];
switch action
    case 'install'
        local_install(app);
    otherwise
        try
            out = local_dispatch(action, app, answers);
        catch ME
            % ea_error already showed a dialog (it leaves the identifier
            % empty); show other errors here, too.
            if ~isempty(ME.identifier)
                local_alert(app, ME.message, 'Permutation testing', 'error');
            end
            ea_cprintf('CmdWinErrors', 'Permutation testing: %s\n', ME.message);
        end
end


%% ------------------------------------------------------------------------
function local_install(app)
if ~isempty(findobj(app.UIFigure, 'Tag', 'ea_permtest_menu'))
    return
end
top = uimenu(app.UIFigure, 'Text', 'Permutation testing', 'Tag', 'ea_permtest_menu');
uimenu(top, 'Text', 'Check settings / preview...', 'MenuSelectedFcn', @(~, ~) ea_permtest_gui_disctract('check', app));
uimenu(top, 'Text', 'Run permutation test...', 'MenuSelectedFcn', @(~, ~) ea_permtest_gui_disctract('run', app));
uimenu(top, 'Text', 'Compute statistics...', 'MenuSelectedFcn', @(~, ~) ea_permtest_gui_disctract('stats', app));
uimenu(top, 'Text', 'Prediction permutation test (k-fold CV)...', 'Separator', 'on', ...
    'MenuSelectedFcn', @(~, ~) ea_permtest_gui_disctract('predict', app));
uimenu(top, 'Text', 'Show significant fibers...', 'Separator', 'on', ...
    'MenuSelectedFcn', @(~, ~) ea_permtest_gui_disctract('show', app));
uimenu(top, 'Text', 'Hide significant fibers', 'MenuSelectedFcn', @(~, ~) ea_permtest_gui_disctract('hide', app));


function out = local_dispatch(action, app, answers)
out = [];
switch action
    case 'check'
        out = local_check(app);
    case 'run'
        out = local_run(app, answers);
    case 'stats'
        out = local_stats(app, answers);
    case 'predict'
        out = local_predict(app, answers);
    case 'show'
        out = local_show(app, answers);
    case 'hide'
        local_hide(app);
    otherwise
        ea_error(['Unknown permutation testing action: ', action]);
end


%% Check settings / preview ------------------------------------------------
function report = local_check(app)
obj = local_obj(app);
d = local_progress(app, 'Checking the current explorer settings...');
cleanup = onCleanup(@() local_close(d));
bundle = ea_permtest_adapter_disctract(obj, local_pathopts(app, obj));
report = ea_permtest_check_disctract(obj, bundle);
clear cleanup

L = {sprintf('Mode: %s', bundle.mode), sprintf('Statistic: %s (%s)', bundle.stattest, bundle.statfamily), ...
    sprintf('Selected patients: %d of %d, mirrored: %d', numel(bundle.patientselection), bundle.Npatients, bundle.mirrorsides), ''};
totalFibers = 0;
for a = 1:numel(bundle.analyses)
    an = bundle.analyses(a);
    fs = bundle.fibersets(an.fiberset);
    nf = cellfun(@numel, fs.nonemptyidx);
    totalFibers = totalFibers + sum(nf);
    nMiss = sum(isnan(bundle.outcomes(an.patientselection, an.outcomeCol)));
    L{end+1} = sprintf('%s: %s fibers, %d missing score(s)', an.label, ...
        strjoin(arrayfun(@(s) sprintf('%s %d', local_side(s), nf(s)), 1:numel(nf), 'Uni', 0), ', '), nMiss); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = sprintf('Disk space of a run with 1000 permutations (R and p): about %.2f GB', totalFibers * 1000 * 4 * 2 / 2^30);
if report.ok
    L{end+1} = 'Consistency check: PASSED (the permutation code reproduces the explorer''s fiber values).';
    local_alert(app, strjoin(L, newline), 'Permutation testing: check', 'success');
else
    L{end+1} = 'Consistency check: FAILED -- see the command window. Do not run permutations with these settings.';
    local_alert(app, strjoin(L, newline), 'Permutation testing: check', 'error');
end


%% Run permutation test ------------------------------------------------------
function run = local_run(app, answers)
obj = local_obj(app);
if isempty(answers)
    name = ['analysis_', datestr(now, 'yyyy-mm-dd')];
    while true
        answers = local_dialog(app, 'Run permutation test', { ...
            'name', 'Analysis name', 'text', name; ...
            'Nperm', 'Number of permutations', 'number', 1000; ...
            'stratifyByGroup', 'Shuffle only within cohorts (patient groups)', 'checkbox', false; ...
            'maxWorkers', 'Parallel workers (only for tests without fast path)', 'dropdown', {'auto', '0 (serial)', '2', '3', '4', '6', '8'}});
        if isempty(answers)
            run = [];
            return
        end
        name = strtrim(answers.name);
        if ~isempty(regexp(name, '^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$', 'once'))
            break
        end
        waitfor(errordlg(sprintf(['Invalid analysis name ''%s''.\n', ...
            'Please use letters, digits, ''-'' and ''_'' only (at most 64 characters).'], name), 'Analysis name', 'modal'));
    end
    answers.name = name;
end
opts = local_pathopts(app, obj);
opts.name = answers.name;
opts.Nperm = answers.Nperm;
opts.stratifyByGroup = logical(answers.stratifyByGroup);
opts.maxWorkers = local_workers(answers.maxWorkers);

d = local_progress(app, sprintf(['Running %d permutations -- progress is shown in the MATLAB command window.\n', ...
    'The explorer is busy until the run has finished.'], opts.Nperm));
cleanup = onCleanup(@() local_close(d));
run = ea_permtest_fibers_disctract(obj, opts);
clear cleanup
local_alert(app, sprintf('Permutation run finished in %.0f s.\n\nResults:\n%s', run.totalSeconds, run.rundir), ...
    'Permutation testing', 'success');


%% Compute statistics -------------------------------------------------------
function out = local_stats(app, answers)
obj = local_obj(app);
out = [];
if isempty(answers)
    [runs, names] = local_runs(app, obj);
    if isempty(runs)
        ea_error('No finished permutation runs found for this explorer -- run a permutation test first.');
    end
    answers = local_dialog(app, 'Compute statistics', { ...
        'rundir', 'Permutation run', 'dropdown', {names, runs}; ...
        'alpha', 'FWE alpha', 'number', 0.05; ...
        'direction', 'Direction', 'dropdown', {{'positive (improvement)', 'negative', 'both, alpha split', 'both, alpha each'}, ...
                                                {'positive', 'negative', 'both-split', 'both-each'}}; ...
        'stepdown', 'Step-down (Westfall-Young)', 'checkbox', true; ...
        'useTailApprox', 'GPD tail approximation', 'checkbox', false});
    if isempty(answers)
        return
    end
end
opts = struct('alpha', answers.alpha, 'direction', answers.direction, 'stepdown', logical(answers.stepdown), ...
    'useTailApprox', logical(answers.useTailApprox));

d = local_progress(app, 'Computing permutation statistics -- see the MATLAB command window.');
cleanup = onCleanup(@() local_close(d));
[stats, files] = ea_permtest_fibers_stats(answers.rundir, opts);
clear cleanup

L = {};
for a = 1:numel(stats.analyses)
    lab = stats.analyses(a).label;
    om = stats.omnibus([stats.omnibus.analysis] == a & strcmp({stats.omnibus.scope}, 'bothsides'));
    ms = stats.maxstat([stats.maxstat.analysis] == a & strcmp({stats.maxstat.family}, 'bothsides'));
    s = sprintf('%s:', lab);
    for o = om
        s = sprintf('%s %s p = %.3g;', s, o.type, o.p);
    end
    for e = ms
        s = sprintf('%s significant fibers (%s, both sides) = %d, step-down %d', s, e.tail, e.NsigSingle, e.NsigStepdown);
    end
    L{end+1} = s; %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = sprintf('%d significant-fiber file(s) exported. Figures and details:', numel(files));
L{end+1} = stats.settings.outdir;
local_alert(app, strjoin(L, newline), 'Permutation statistics', 'success');
if local_hasui(app)
    setappdata(app.UIFigure, 'ea_permtest_lastexports', strrep(stats.settings.outdir, ...
        [filesep, 'stats', filesep], [filesep, 'exports', filesep]));
end
out = struct('stats', stats, 'files', {files});


%% Prediction permutation test (k-fold CV) ---------------------------------
function run = local_predict(app, answers)
obj = local_obj(app);
run = [];
[kIter, posthoc] = local_cvpanel(app, obj);
if isempty(answers)
    [~, names] = local_runs(app, obj);
    names = names(~startsWith(names, '(older run)'));
    if isempty(names)
        name = ['analysis_', datestr(now, 'yyyy-mm-dd')];
    else
        name = names{1}; % newest analysis of this explorer: add the prediction run to it
    end
    targets = [{'(none)'}, obj.M.clinical.labels(:)'];
    fields = { ...
        'info', sprintf('From the CV panel: k = %d, iterations = %d, post-hoc correction for group: %s', ...
            obj.kfold, kIter, char(string(posthoc))), 'label', []; ...
        'name', 'Analysis name', 'text', name; ...
        'Nperm', 'Number of permutations', 'number', 1000; ...
        'stratifyByGroup', 'Shuffle only within cohorts (patient groups)', 'checkbox', false};
    if strcmp(obj.multitractmode, 'Split & Color By Subscore')
        fields(end+1, :) = {'targetLabel', 'Combined prediction vs (optional)', 'dropdown', targets};
    end
    while true
        answers = local_dialog(app, 'Prediction permutation test (k-fold CV)', fields);
        if isempty(answers)
            return
        end
        name = strtrim(answers.name);
        if ~isempty(regexp(name, '^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$', 'once'))
            break
        end
        fields{2, 4} = name;
        waitfor(errordlg(sprintf(['Invalid analysis name ''%s''.\n', ...
            'Please use letters, digits, ''-'' and ''_'' only (at most 64 characters).'], name), 'Analysis name', 'modal'));
    end
    answers.name = name;
    if ~isfield(answers, 'targetLabel') || strcmp(answers.targetLabel, '(none)')
        answers.targetLabel = '';
    end
    answers.kIter = kIter;
    answers.posthoccorrect = posthoc;
end
opts = local_pathopts(app, obj);
for f = {'name', 'Nperm', 'targetLabel', 'kIter', 'posthoccorrect'}
    if isfield(answers, f{1})
        opts.(f{1}) = answers.(f{1});
    end
end
opts.stratifyByGroup = logical(answers.stratifyByGroup);
if local_hasui(app)
    opts.confirmFcn = @(msg) strcmp(questdlg(msg, 'Prediction permutation test', 'Continue', 'Cancel', 'Cancel'), 'Continue');
end

d = local_progress(app, ['Prediction permutation test -- progress is shown in the MATLAB command window.', newline, ...
    'The explorer is busy until the run has finished.']);
cleanup = onCleanup(@() local_close(d));
run = ea_permtest_cv_disctract(obj, opts);
clear cleanup
if run.cancelled
    return
end
r = run.results;
L = {};
for v = 1:numel(r.voterLabels)
    L{end+1} = sprintf('%s: Spearman R = %.3f (p = %.3g), Pearson R = %.3f (p = %.3g)', ...
        r.voterLabels{v}, r.Rvoter(v, 1), r.pvoter(v, 1), r.Rvoter(v, 2), r.pvoter(v, 2)); %#ok<AGROW>
end
if isfield(r, 'Rcomb')
    L{end+1} = sprintf('Combined vs %s: Spearman R = %.3f (p = %.3g), Pearson R = %.3f (p = %.3g)', ...
        r.targetLabel, r.Rcomb(1), r.pcomb(1), r.Rcomb(2), r.pcomb(2));
end
L{end+1} = '';
L{end+1} = run.rundir;
local_alert(app, strjoin(L, newline), 'Prediction permutation test', 'success');


function [kIter, posthoc] = local_cvpanel(app, obj)
% k-fold iterations and post-hoc correction as set in the explorer's CV panel
kIter = obj.kIter;
posthoc = false;
if isobject(app) && isprop(app, 'IterationsEditField') && isvalid(app.IterationsEditField)
    kIter = app.IterationsEditField.Value;
end
if isobject(app) && isprop(app, 'PosthoccorrectforgroupCheckBox') && isvalid(app.PosthoccorrectforgroupCheckBox)
    posthoc = logical(app.PosthoccorrectforgroupCheckBox.Value);
end


%% Show / hide significant fibers ------------------------------------------
function h = local_show(app, answers)
obj = local_obj(app);
if isempty(answers)
    start = local_lastexports(app, obj);
    [f, p] = uigetfile('*.mat', 'Choose a significant-fiber export file', start);
    if isequal(f, 0)
        h = [];
        return
    end
    answers = struct('file', fullfile(p, f));
end
E = load(answers.file);
if ~isfield(E, 'fibcell') || ~isfield(E, 'vals')
    ea_error('This file is not a significant-fiber export (fibcell / vals missing).');
end
if isempty(obj.resultfig) || ~isgraphics(obj.resultfig)
    ea_error('The 3D figure of the explorer is closed -- recreate it (Tools > Recreate visualization figure) first.');
end
local_hide(app);

allv = vertcat(E.vals{:});
vmax = max(abs(allv));
prefs = ea_prefs;
set(0, 'CurrentFigure', obj.resultfig);
h = {};
for side = 1:numel(E.fibcell)
    if isempty(E.fibcell{side})
        continue
    end
    t = streamtube(E.fibcell{side}, prefs.d3.fiberwidth);
    nones = repmat({'none'}, size(t));
    [t.EdgeColor] = nones{:};
    cols = arrayfun(@(v) local_color(v, vmax), E.vals{side}(:), 'Uni', 0);
    [t.FaceColor] = cols{:};
    h{end+1} = t; %#ok<AGROW>
end
setappdata(obj.resultfig, 'ea_permtest_sigfibers', h);
n = sum(cellfun(@numel, E.fibcell));
fprintf('Showing %d significant fibers from %s\n', n, answers.file);
if isfield(E, 'info')
    fprintf('  %s, %s, %s, %s, alpha %g\n', E.info.analysis, E.info.family, E.info.direction, E.info.method, E.info.alpha);
end


function local_hide(app)
obj = local_obj(app);
if isempty(obj.resultfig) || ~isgraphics(obj.resultfig)
    return
end
h = getappdata(obj.resultfig, 'ea_permtest_sigfibers');
for i = 1:numel(h)
    delete(h{i}(isgraphics(h{i})));
end
setappdata(obj.resultfig, 'ea_permtest_sigfibers', {});


function c = local_color(v, vmax)
% Light to dark red for positive R, light to dark blue for negative R.
w = min(abs(v) / max(vmax, eps), 1);
if v >= 0
    c = (1 - w) * [1.00 0.75 0.60] + w * [0.70 0.05 0.05];
else
    c = (1 - w) * [0.65 0.80 1.00] + w * [0.05 0.20 0.70];
end


%% Helpers ------------------------------------------------------------------
function obj = local_obj(app)
obj = app.tractset;
if isempty(obj) || ~isa(obj, 'ea_disctract')
    ea_error('No fiber filtering analysis is loaded in the explorer.');
end


function pathopts = local_pathopts(app, obj)
% Lead group folder for the shuffle set and the results.
if ~isempty(obj.analysispath) && isfile(obj.analysispath)
    pathopts = struct('fibfiltPath', obj.analysispath);
else
    pathopts = struct('groupdir', local_groupdir(app, obj));
end


function groupdir = local_groupdir(app, obj)
% The explorer keeps analysispath up to date (it is set when a .fibfilt is
% opened or saved): <groupdir>/fiberfiltering/<ID>.fibfilt. Without it,
% the lead group path stored in the explorer can be outdated (e.g. another
% machine), so the user is asked once per explorer session.
if ~isempty(obj.analysispath) && isfile(obj.analysispath)
    groupdir = fileparts(fileparts(obj.analysispath));
    return
end
key = 'ea_permtest_groupdir';
if local_hasui(app) && isappdata(app.UIFigure, key)
    stored = getappdata(app.UIFigure, key);
    if isfield(stored, 'ID') && strcmp(stored.ID, obj.ID) && isfolder(stored.groupdir)
        groupdir = stored.groupdir;
        return
    end
end
groupdir = fileparts(obj.leadgroup);
if ~local_hasui(app)
    return % command line: the adapter checks the folder
end
start = groupdir;
if ~isfolder(start)
    start = pwd;
end
chosen = uigetdir(start, 'Lead group folder (permutation results are saved below it)');
if isequal(chosen, 0)
    ea_error('No lead group folder chosen -- nothing was calculated.', 'title', 'Permutation testing');
end
groupdir = chosen;
setappdata(app.UIFigure, key, struct('ID', obj.ID, 'groupdir', groupdir));


function root = local_resultsroot(app, obj)
root = fullfile(local_groupdir(app, obj), 'fiberfiltering', 'permutation');


function [runs, names] = local_runs(app, obj)
% Finished fiber runs of this explorer, newest first:
%   <root>/<analysis name>/fibers         (named analyses of this explorer ID)
%   <root>/<ID>/<ID>_fiberperm_*          (runs made before named analyses)
root = local_resultsroot(app, obj);
runs = {};
names = {};
dates = [];
D = dir(fullfile(root, '*', 'fibers', 'run_info.mat'));
for i = 1:numel(D)
    f = D(i).folder;
    if ~ismember('seconds', who('-file', fullfile(f, 'run_info.mat')))
        continue % not finished
    end
    S = load(fullfile(f, 'run_info.mat'), 'settings');
    if ~isfield(S.settings, 'ID') || ~strcmp(S.settings.ID, obj.ID)
        continue
    end
    [~, name] = fileparts(fileparts(f));
    runs{end+1} = f; %#ok<AGROW>
    names{end+1} = name; %#ok<AGROW>
    dates(end+1) = D(i).datenum; %#ok<AGROW>
end
D = dir(fullfile(root, obj.ID, '*_fiberperm_*'));
D = D([D.isdir]);
for i = 1:numel(D)
    f = fullfile(D(i).folder, D(i).name);
    if isfile(fullfile(f, 'run_info.mat')) && ismember('seconds', who('-file', fullfile(f, 'run_info.mat')))
        runs{end+1} = f; %#ok<AGROW>
        names{end+1} = ['(older run) ', regexprep(D(i).name, ['^', regexptranslate('escape', obj.ID), '_'], '')]; %#ok<AGROW>
        dates(end+1) = D(i).datenum; %#ok<AGROW>
    end
end
[~, o] = sort(dates, 'descend');
runs = runs(o);
names = names(o);


function start = local_lastexports(app, obj)
start = '';
if local_hasui(app) && isappdata(app.UIFigure, 'ea_permtest_lastexports')
    start = getappdata(app.UIFigure, 'ea_permtest_lastexports');
end
if isempty(start) || ~isfolder(start)
    runs = local_runs(app, obj);
    if ~isempty(runs) && isfolder(fullfile(runs{1}, 'exports'))
        start = fullfile(runs{1}, 'exports');
    elseif ~isempty(runs)
        start = runs{1};
    else
        start = local_resultsroot(app, obj);
    end
end
if isfolder(start)
    start = [start, filesep];
end


function w = local_workers(v)
if isnumeric(v)
    w = v;
elseif strcmp(v, 'auto')
    w = 'auto';
else
    w = str2double(strtok(v));
end


function s = local_side(side)
labels = {'right', 'left'};
s = labels{min(side, 2)};


function ok = local_hasui(app)
% true when the explorer window is available for dialogs
ok = (isobject(app) && isprop(app, 'UIFigure')) || (isstruct(app) && isfield(app, 'UIFigure'));
ok = ok && ~isempty(app.UIFigure) && isgraphics(app.UIFigure);


function local_alert(app, msg, title, icon)
shown = false;
if local_hasui(app)
    try
        uialert(app.UIFigure, msg, title, 'Icon', icon);
        shown = true;
    catch
    end
end
if ~shown
    fprintf('[%s] %s\n', title, msg);
end


function d = local_progress(app, msg)
d = [];
if local_hasui(app)
    try
        d = uiprogressdlg(app.UIFigure, 'Title', 'Permutation testing', 'Message', msg, 'Indeterminate', 'on');
        drawnow;
    catch
        d = [];
    end
end


function local_close(d)
if ~isempty(d) && isvalid(d)
    close(d);
end


function answers = local_dialog(app, titlestr, fields)
% Small modal dialog. fields: rows of {name, label, type, default}, type
% 'number' | 'text' | 'checkbox' | 'dropdown' (default = {items} or {items, values})
% | 'label' (information only).
answers = [];
n = size(fields, 1);
rowH = 32;
fig = uifigure('Name', titlestr, 'Position', [0 0 520, rowH * n + 70], 'WindowStyle', 'modal', 'Resize', 'off');
movegui(fig, 'center');
if local_hasui(app)
    p = app.UIFigure.Position;
    fig.Position(1:2) = p(1:2) + (p(3:4) - fig.Position(3:4)) / 2;
end
ctrl = cell(1, n);
for i = 1:n
    y = fig.Position(4) - 20 - i * rowH;
    switch fields{i, 3}
        case 'label'
            ctrl{i} = uilabel(fig, 'Text', fields{i, 2}, 'Position', [20 y 480 22], 'WordWrap', 'on');
        case 'checkbox'
            ctrl{i} = uicheckbox(fig, 'Text', fields{i, 2}, 'Value', fields{i, 4}, 'Position', [20 y 480 22]);
        otherwise
            uilabel(fig, 'Text', fields{i, 2}, 'Position', [20 y 250 22]);
            if strcmp(fields{i, 3}, 'number')
                ctrl{i} = uieditfield(fig, 'numeric', 'Value', fields{i, 4}, 'Position', [280 y 220 22]);
            elseif strcmp(fields{i, 3}, 'text')
                ctrl{i} = uieditfield(fig, 'text', 'Value', fields{i, 4}, 'Position', [280 y 220 22]);
            else
                def = fields{i, 4};
                if iscell(def{1})
                    ctrl{i} = uidropdown(fig, 'Items', def{1}, 'ItemsData', def{2}, 'Position', [280 y 220 22]);
                else
                    ctrl{i} = uidropdown(fig, 'Items', def, 'Position', [280 y 220 22]);
                end
            end
    end
end
uibutton(fig, 'Text', 'Cancel', 'Position', [300 15 95 26], 'ButtonPushedFcn', @(~, ~) uiresume(fig));
uibutton(fig, 'Text', 'OK', 'Position', [405 15 95 26], 'ButtonPushedFcn', @(~, ~) local_ok(fig));
uiwait(fig);
if isvalid(fig)
    ok = isappdata(fig, 'ok');
    if ok
        answers = struct;
        for i = 1:n
            if ~strcmp(fields{i, 3}, 'label')
                answers.(fields{i, 1}) = ctrl{i}.Value;
            end
        end
    end
    delete(fig);
end


function local_ok(fig)
setappdata(fig, 'ok', true);
uiresume(fig);
