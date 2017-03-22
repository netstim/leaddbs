function options=ea_getptopts(directory,options)
% function to generate minimal options struct from a patient directory.
% prior options struct can (optionally) be fed in and only patient specific
% entries will be replaced.
if isempty(directory)
    directory=pwd;
end

if strcmp(directory(end),filesep) % strip trailing filesep
    directory=directory(1:end-1);
end

[options.root,options.patientname]=fileparts(directory);
options.root=[options.root,filesep];

options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
