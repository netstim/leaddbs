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
    for fi=f
        % apply biasfieldcorrection and reorient/crop
        ea_anatpreprocess([directory,options.prefs.rawpreniis{fi}]);
    end
    try
    fs=fopen([directory,'.pp'],'w');
    fprintf(fs,'%s','dcmbiasfielddone');
    fclose(fs);
    end
end

try
    options.prefs.prenii_unnormalized=options.prefs.rawpreniis{f(1)};
    options.primarytemplate=ptemps{f(1)};
catch
    options.primarytemplate='';
end

presentfiles=options.prefs.rawpreniis(f);


