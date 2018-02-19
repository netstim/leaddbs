function home = ea_gethome

if ispc
    home = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    home = getenv('HOME');
end

home = [home, filesep];

% TP: Bug when home is pointing to a drive: i.e. C:\ then another filesep
% is added at end to make it C:\\ causing potential error. Simple fix below:
if strcmp(home(end-1), filesep)
    home(end) = [];
end

if isdeployed
   mkdir([ctfroot, filesep, 'home', filesep]);
   home=[ctfroot, filesep, 'home', filesep];
end
