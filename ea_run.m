function ea_run(cmd,options)
% This function is the main execution function of Lead-DBS. It is
% distributed within Lead-DBS toolbox (www.lead-dbs.org)
% __________________________________________________________________________________
% Copyright (C) 2016 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ~isfield(options,'lc') % might be predefined from an exported script..
    try
        options.lc=load([fileparts(which('lead')),filesep,'connectomics',filesep,'lc_options.mat']);
    catch
        options.lc=[];
    end
end

if options.d3.autoserver && options.d3.write
    choice = questdlg('Are you sure you want to export results to the server automatically?', ...
        'Auto Server-Export', ...
        'Cancel','Yes','Cancel');
    % Handle response
    switch choice
        case 'Cancel'
            return
    end
end

clc
uipatdirs=options.uipatdirs;

if isempty(uipatdirs)
    uipatdirs={'No Patient Selected'};
end

prefs=ea_prefs('');
if length(uipatdirs)>1 && ~isempty(which('parpool')) && prefs.pp.do && ~strcmp(cmd,'export') % do parallel processing if available and set in ea_prefs.
    try delete(gcp); end
    pp=parpool(prefs.pp.profile,prefs.pp.csize);
    
    for pat=1:length(uipatdirs)
        % set patient specific options
        opts{pat}=options;
        opts{pat}.root=[fileparts(uipatdirs{pat}),filesep];
        [~,thispatdir]=fileparts(uipatdirs{pat});
        opts{pat}.patientname=thispatdir;
    end
    
    parfor pat=1:length(uipatdirs)
        
        % run main function
        try
            switch cmd
                case 'run'
                    ea_autocoord(opts{pat});
            end
        catch
            warning([opts{pat}.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.' ]);
        end
        
    end
    delete(pp);
    
else
    switch cmd
        
        case 'export'
            
            ea_export(options);
            
        case 'run'
            
            for pat=1:length(uipatdirs)
                % set patient specific options
                options.root=[fileparts(uipatdirs{pat}),filesep];
                [root,thispatdir]=fileparts(uipatdirs{pat});
                options.patientname=thispatdir;
                % run main function
                
                if length(uipatdirs)>1 % multi mode. Dont stop at errors.
                    try
                                ea_autocoord(options);
                    catch
                        warning([options.patientname,' failed. Please run this patient again and adjust parameters. Moving on to next patient.' ]);
                    end
                else
                            ea_autocoord(options);
                end
            end
            
    end
end





