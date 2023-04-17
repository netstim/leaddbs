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

cmd = ea_lcm_resolvecmd(options.lcm.cmd);
[sfile]=ea_handleseeds(options.lcm.seeds);
if strcmp(cmd,'seed')
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

    if regexp(options.lcm.func.connectome, '^Patient''s fMRI - ', 'once')
        % native space nifti file
        options.root=[options.uivatdirs{1}];
        [options.root,options.patientname]=fileparts(options.root);
        options.root=[options.root,filesep];
        restfname = options.lcm.func.connectome(length('Patient''s fMRI - ')+1:end);
        seed_tc=cs_fmri_conseed_nifti([options.uivatdirs{1},filesep,restfname,'.nii'],tsfile,options);
        [pth,fname,ext]=fileparts(tsfile{1});
        save(fullfile(pth,[fname,'.mat']),'seed_tc');
    else
        cname=options.lcm.func.connectome;
        if ismember('>',cname)
            delim=strfind(cname,'>');
            subset=cname(delim+1:end);
            cname=cname(1:delim-1);
        end
        dataset=loadjson([ea_getconnectomebase,'fMRI',filesep,cname,filesep,'dataset_info.json']);
        switch cmd
            case 'seed'
                switch dataset.type
                    case 'fMRI_matrix'
                        cs_fmri_conseed_seed_matrix(ea_getconnectomebase,...
                            options.lcm.func.connectome,...
                            tsfile,...
                            cmd,...
                            '0',...
                            options.lcm.odir,...
                            options.lcm.omask);
                    case 'fMRI_timecourses'
                        cs_fmri_conseed_seed_tc(ea_getconnectomebase,...
                            options.lcm.func.connectome,...
                            tsfile,...
                            cmd,...
                            '0',...
                            options.lcm.odir,...
                            options.lcm.omask);
                end
            case 'pseed'
                 switch dataset.type
                    case 'fMRI_matrix'
                        ea_error('Command partial seed is not supported for matrix type datasets.')
                    case 'fMRI_timecourses'
                        cs_fmri_conseed_pseed(ea_getconnectomebase,...
                            options.lcm.func.connectome,...
                            tsfile,...
                            cmd,...
                            '0',...
                            options.lcm.odir,...
                            options.lcm.omask);
                 end
            case {'matrix','pmatrix'}
                switch dataset.type
                    case 'fMRI_matrix'
                        cs_fmri_conseed_matrix_matrix(ea_getconnectomebase,...
                            options.lcm.func.connectome, ...
                            tsfile,...
                            options.lcm.odir);
                    case 'fMRI_timecourses'
                        cs_fmri_conseed_matrix(ea_getconnectomebase,...
                            options.lcm.func.connectome,...
                            tsfile,...
                            cmd,...
                            '0',...
                            options.lcm.odir,...
                            options.lcm.func.exportgmtc);
                end
            case 'pmap'
                switch dataset.type
                    case 'fMRI_matrix'
                        ea_error('Command partial map is not supported for matrix type datasets.');
                    case 'fMRI_timecourses'
                        cs_fmri_conseed_pmap(ea_getconnectomebase,...
                            options.lcm.func.connectome,...
                            tsfile,...
                            cmd,...
                            '0',...
                            options.lcm.odir,...
                            options.lcm.omask);
                end
        end
    end
end
disp('Done.');
