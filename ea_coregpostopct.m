function ea_coregpostopct(options)
% Entry function to coregister post-op CT to pre-op MRI

if ~ea_reglocked(options, options.subj.postopAnat.CT.coreg)
    % Setup log
    ea_mkdir(fileparts(options.subj.coreg.log.logBaseName));
    diary([options.subj.coreg.log.logBaseName, 'CT', datestr(now, 'yyyymmddTHHMMss'), '.log']);

    % Do coregistration
    switch lower(options.coregct.method)
        case lower({'ANTs (Avants 2008)', 'ANTs'})
            ea_coregpostopct_ants(options);
        case lower({'ANTs (Avants 2008) multiple runs', 'ANTsMulti'})
            ea_coregpostopct_ants_multiple(options);
        case lower({'ANTs (Avants 2008) multiple runs + subcortical refine', 'ANTsMultiScrf'})
            ea_coregpostopct_ants_multiple_refine(options);
        case lower({'ANTs (Avants 2008) + subcortical refine', 'ANTsScrf'})
            ea_coregpostopct_ants_refine(options);
        case lower({'BRAINSFit (Johnson 2007)', 'BRAINSFit'})
            ea_coregpostopct_brainsfit(options);
        case lower({'FLIRT (Jenkinson 2001 & 2002)', 'FLIRT'})
            ea_coregpostopct_fsl(options);
        otherwise
            warning('Coregistrion method not recognized...');
            diary off;
            return;
    end

    % Dump method
    ea_dumpmethod(options, 'coreg');

    % Compute tone-mapped coregistered CT
    if options.modality == 2
        ea_tonemapct_file(options, 'native');
    end

    diary off;
end
