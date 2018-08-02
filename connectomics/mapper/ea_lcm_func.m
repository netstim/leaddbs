function ea_lcm_func(options)

if strcmp(options.lcm.func.connectome,'No functional connectome found.')
    return
end

disp('Running functional connectivity...');

if iscell(options.lcm.seeds)
   if length(options.lcm.seeds)==1
      options.lcm.seeds=options.lcm.seeds{1}; % supply as char
   end
end

[sfile]=ea_handleseeds(options.lcm.seeds);
if strcmp(ea_lcm_resolvecmd(options.lcm.cmd),'seed')
    chunk=options.prefs.lcm.chunk;
    if ~chunk % split up seedlist for memory considerations.
        chunk=length(sfile);
    end
else
    chunk=length(sfile); % need to consider all seeds in one go for all other commands except seed.
end

for run=1:chunk:length(sfile)
    try
        tsfile=sfile(run:run+chunk-1);
    catch
        tsfile=sfile(run:end);
    end

    if isfield(options,'uivatdirs')
        if ~isempty(options.uivatdirs)
            options.lcm.odir=[];
        end
    end

    if strncmp(options.lcm.func.connectome, 'Patient-specific fMRI - ', 24)
        % native space nifti file
        options.root=[options.uivatdirs{1}];
        [options.root,options.patientname]=fileparts(options.root);
        options.root=[options.root,filesep];
        restfname = options.lcm.func.connectome(25:end);
        cs_fmri_conseed_nifti([options.uivatdirs{1},filesep,restfname,'.nii'],tsfile,options)
    else
        cs_fmri_conseed(ea_getconnectomebase,options.lcm.func.connectome,...
            tsfile,...
            ea_lcm_resolvecmd(options.lcm.cmd),...
            '0',...
            options.lcm.odir,...
            options.lcm.omask);
    end
end
disp('Done.');
