function options = ea_get_OSS_DBS_options(patientPath, N_tracts, ConnectomeName, native)

    % update options for OSS-DBS launch without GUI
    % unupdated parameters are the same as during the last GUI settings session

    options = ea_prefs('');
    options = ea_getptopts(patientPath, options);
   
    options.native = native;
    options.groupmode = 1;
    options.groupid = 'currentune';
    options.stimSetMode = 0;  % we do not use it for current-optimization
    options.orignative = 1;
    options.prefs.machine.vatsettings.butenko_interactive = 0;
    options.prefs.machine.vatsettings.butenko_calcAxonActivation = 1;
    %options.prefs.machine.vatsettings.butenko_connectome = ['Multi-Tract: ', ConnectomeName];

    % default settings for axons 
    %options.prefs.machine.vatsettings.butenko_axonLength = 15.0 * ones(N_tracts,1);   % 15 mm, smaller axons will not be seeded!
    %options.prefs.machine.vatsettings.butenko_fiberDiameter = 3.0 * ones(N_tracts,1);   % 3 micro m 

    options.earoot=ea_getearoot;
    options.verbose=3;
    options.sides= 1;  % does this change anything?
    options.fiberthresh=1;
    options.writeoutstats=1;
    options.writeoutpm = 0;
    options.colormap='parula(64)';
    options.d3.write=1;
    options.d3.prolong_electrode=2;
    options.d3.writeatlases=1;
    options.macaquemodus=0;    
end
