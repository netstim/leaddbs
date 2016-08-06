function [options,presentfiles]=ea_assignpretra(options)


directory=[options.root,options.patientname,filesep];
ptemps={'','_t1','_pd'};
for prenii=1:length(options.prefs.rawpreniis)

    fexist(prenii)=exist([directory,options.prefs.rawpreniis{prenii}],'file');
    
end

f=find(fexist);
if isempty(f)
    warning(['No anatomy information found. Please put either ',options.prefs.rawpreniis{1},', ',  options.prefs.rawpreniis{2},' or ', options.prefs.rawpreniis{3}, ' into subject folder.']);
end


if ~exist([directory,'.pp'],'file') && ~exist([directory,'ea_normmethod_applied.mat'],'file') % only do this once, small hidden flag .pp inside patient folder will show this has been done before.
    for f=find(fexist)
        % apply biasfieldcorrection and reorient/crop
        ea_dcm2nii([directory,options.prefs.rawpreniis{f}]);
        ea_bias_field_correction([directory,options.prefs.rawpreniis{f}])
    end
    fs=fopen([directory,'.pp'],'w');
end
try
    options.prefs.prenii_unnormalized=options.prefs.rawpreniis{f(1)};
    options.primarytemplate=ptemps{f(1)};
catch
    options.primarytemplate=options.prefs.prenii_unnormalized;
end

presentfiles=options.prefs.rawpreniis(f);


