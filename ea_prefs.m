function prefs=ea_prefs(patientname)
%#ea_prefs_default
if ~exist('patientname','var')
    patientname='';
end

% get default prefs
prefs=ea_prefs_default(patientname);
% now overwrite with user prefs stored in /home
home=ea_gethome;
uid=['ea_prefs_',dash2sub(ea_generate_guid)];
if isdeployed
    disp(['Running Lead-DBS in compiled mode, CTFROOT=',ea_getearoot,'; HOME=',home,'.']);
end

if ~exist([home,'.ea_prefs.m'],'file')
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.m'],[home,'.ea_prefs.m']);
end
if ~exist([home,'.ea_prefs.mat'],'file')
    copyfile([ea_getearoot,'common',filesep,'ea_prefs_default.mat'],[home,'.ea_prefs.mat']);
end


defmachine=load([ea_getearoot,'common',filesep,'ea_prefs_default.mat']);
defmachine=defmachine.machine;
try
    copyfile([home,'.ea_prefs.m'],[ea_getearoot,uid,'.m'])
    uprefs=feval(uid,patientname);
    delete([ea_getearoot,uid,'.m']);
    load([home,'.ea_prefs.mat']);
catch
   warning('User preferences file could not be read. Please set write permissions to Lead-DBS install directory accordingly.');
   return
end

prefs=combinestructs(prefs,uprefs);

machine=combinestructs(defmachine,machine);

prefs.machine=machine;



function str=dash2sub(str) % replaces subscores with spaces
str(str=='-')='_';


function prefs=combinestructs(prefs,uprefs)

ufn = fieldnames(uprefs);

% Handle struct arrays of different sizes
if length(uprefs) > length(prefs)
    prefs = cat(prefs, uprefs(length(prefs)+1:end));
end

% Iterate through struct array. Most prefs are a single struct, not an array.
for sa_ix = 1:length(uprefs)
    for fn_ix = 1:length(ufn)
        fn = ufn{fn_ix};
        if isstruct(uprefs(sa_ix).(fn)) && isfield(prefs(sa_ix),fn)
            prefs(sa_ix).(fn) = combinestructs(prefs(sa_ix).(fn), uprefs(sa_ix).(fn));
        else
            prefs(sa_ix).(fn) = uprefs(sa_ix).(fn);
        end
    end
end
