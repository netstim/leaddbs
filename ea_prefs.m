function prefs = ea_prefs(patientname, context)

if ~exist('patientname', 'var')
    patientname = '';
end

% 'normal' context by default. If context is not 'normal' (cluster run or
% exported job case), 'machine' settings will not be added to 'prefs' here.
% The value will be read from options.prefs.machines.
if ~exist('context', 'var')
    context = 'normal';
end

% get default prefs
dprefs = ea_prefs_default(patientname);
dmachine = load([ea_getearoot, 'common', filesep, 'ea_prefs_default.mat']);

% Check default auto colormap setting
try
    eval([dprefs.d3.roi.defaultcolormap, ';'])
catch
    ea_cprintf('CmdWinWarnings', 'Default auto colormap was not set properly, fallback to ''turbo'' now.\n');
    ea_setprefs('d3.roi.defaultcolormap', 'turbo', 'user');
end

% if isdeployed
%     disp(['Running Lead-DBS in compiled mode, CTFROOT=', ea_getearoot, '; HOME=', ea_gethome, '.']);
% end

if ~isfile(ea_prefspath(ea_prefsext))
    copyfile([ea_getearoot, 'common', filesep, 'ea_prefs_default', ea_prefsext], ea_prefspath(ea_prefsext), 'f');
end

if ~isfile(ea_prefspath('mat'))
    copyfile([ea_getearoot, 'common', filesep, 'ea_prefs_default.mat'], ea_prefspath('mat'), 'f');
end

% load user prefs
try
    if ~isdeployed
        % Temporarily add ~/.leaddbs to search path and get user prefs.
        addpath(ea_prefsdir);
        uprefs = ea_prefs_user(patientname);
        rmpath(ea_prefsdir);
        umachine = load(ea_prefspath('mat'));
    else
        fid = fopen(ea_prefspath('json'),'rt');
        uprefs = jsondecode(fread(fid,'*char')'); fclose(fid);
        umachine = load(ea_prefspath('mat'));
    end
catch ME
    prefs=dprefs; % it seems user-defined prefs cannot be loaded.
    warning(ME.message);
    return
end

% combine default prefs and user prefs
prefs = combinestructs(dprefs, uprefs);

% add machine prefs only when lead task is not running on a cluster or with
% exported job file. In both cases, machine prefs will be defined in
% options.prefs.machine.
if strcmp(context, 'normal')
    prefs.machine = combinestructs(dmachine.machine, umachine.machine);
end

% legacy code support for gl/l normalized file differentiation:
prefs.prenii=prefs.gprenii;
prefs.tranii=prefs.gtranii;
prefs.cornii=prefs.gcornii;
prefs.sagnii=prefs.gsagnii;
prefs.ctnii=prefs.gctnii;


function prefs = combinestructs(prefs, uprefs)

ufn = fieldnames(uprefs);

% Handle struct arrays of different sizes
if length(uprefs) > length(prefs)
    prefs = cat(prefs, uprefs(length(prefs)+1:end));
end

% Iterate through struct array. Most prefs are a single struct, not an array.
for sa_ix = 1:length(uprefs)
    for fn_ix = 1:length(ufn)
        fn = ufn{fn_ix};
        if isstruct(uprefs(sa_ix).(fn)) && isfield(prefs(sa_ix), fn)
            prefs(sa_ix).(fn) = combinestructs(prefs(sa_ix).(fn), uprefs(sa_ix).(fn));
        else
            prefs(sa_ix).(fn) = uprefs(sa_ix).(fn);
        end
    end
end
