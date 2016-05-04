function options=ea_assignpretra(options)


directory=[options.root,options.patientname,filesep];
ptemps={'','_t1','_pd'};
for prenii=1:length(options.prefs.rawpreniis)

    fexist(prenii)=exist([directory,options.prefs.rawpreniis{prenii}],'file');
    
end

f=find(fexist);
if isempty(f)
    ea_error(['No anatomy information found. Please put either ',options.prefs.rawpreniis{1},', ',  options.prefs.rawpreniis{2},' or ', options.prefs.rawpreniis{3}, ' into subject folder.']);
end
options.prefs.prenii_unnormalized=options.prefs.rawpreniis{f(1)};
options.primarytemplate=ptemps{f(1)};


