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


ufn=fieldnames(uprefs);

for uf=1:length(ufn) % compare user preferences with defaults and overwrite defaults where present.
    if isstruct(uprefs.(ufn{uf}))
        ufn2=fieldnames(uprefs.(ufn{uf}));
        for uf2=1:length(ufn2)
            if isstruct(uprefs.(ufn{uf}).(ufn2{uf2}))
                ufn3=fieldnames(uprefs.(ufn{uf}).(ufn2{uf2}));
                for uf3=1:length(ufn3)
                    if isstruct(uprefs.(ufn{uf}).(ufn2{uf2}).(ufn3{uf3}));
                        ufn4=fieldnames(uprefs.(ufn{uf}).(ufn2{uf2}).(ufn3{uf3}));
                        for uf4=1:length(ufn4) % add fourth level entries
                            prefs.(ufn{uf}).(ufn2{uf2}).(ufn3{uf3}).(ufn4{uf4})=uprefs.(ufn{uf}).(ufn2{uf2}).(ufn3{uf3}).(ufn4{uf4});
                        end
                    else % add third level entries
                        prefs.(ufn{uf}).(ufn2{uf2}).(ufn3{uf3})=uprefs.(ufn{uf}).(ufn2{uf2}).(ufn3{uf3});
                    end
                end
            else % add second level entries
                prefs.(ufn{uf}).(ufn2{uf2})=uprefs.(ufn{uf}).(ufn2{uf2});
            end
        end
    else % add first level entries
        prefs.(ufn{uf})=uprefs.(ufn{uf});
    end
end