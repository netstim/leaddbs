function prefs=ea_prefs(patientname)


% get default prefs
prefs=ea_prefs_default(patientname);
% now overwrite with user prefs stored in /home
home=ea_gethome;
uid=['ea_prefs_',dash2sub(ea_generate_guid)];

if ~exist([home,'.ea_prefs.m'],'file');
    copyfile([ea_getearoot,'ea_prefs_default.m'],[home,'.ea_prefs.m']);
end
try
    copyfile([home,'.ea_prefs.m'],[ea_getearoot,uid,'.m'])
    uprefs=feval(uid,patientname);
    delete([ea_getearoot,uid,'.m']);
end

ufn=fieldnames(uprefs);

for uf=1:length(ufn)
   prefs.(ufn{uf})=uprefs.(ufn{uf}); 
end



function str=dash2sub(str) % replaces subscores with spaces
str(str=='-')='_';
