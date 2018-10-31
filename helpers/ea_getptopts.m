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
options.modality=ea_getmodality(directory);
[options.root,options.patientname]=fileparts(directory);
if isempty(options.root)
    options.root=[pwd,filesep];
else
    options.root=[options.root,filesep];
end

options.prefs=ea_prefs(options.patientname);
options=ea_assignpretra(options);
options.earoot=ea_getearoot;
options=ea_detsides(options);
if exist([directory,filesep,'ea_reconstruction.mat'],'file')
    load([directory,filesep,'ea_reconstruction.mat']);
    options.elmodel=reco.props(1).elmodel;
end

