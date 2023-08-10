function ea_normalize(options)
% Entry function to run normalization

if ~ea_reglocked(options, options.subj.preopAnat.(options.subj.AnchorModality).norm)
    % Setup log
    if options.prefs.diary
        ea_mkdir(fileparts(options.subj.norm.log.logBaseName));
        diary([options.subj.norm.log.logBaseName, datestr(now, 'yyyymmddTHHMMss'), '.log']);
    end

    % Dump method
    % Do it before runnung normalization, needed by applying norm functions
    if ~ismember(lower(options.normalize.method), lower({'(Re-)apply (priorly) estimated normalization', 'Apply', 'ReApply'}))
        ea_dumpmethod(options, 'norm');
    end

    % Do normalization
    switch lower(options.normalize.method)
        case lower({'ANTs (Avants 2008)', 'ANTs'})
            ea_normalize_ants(options);
        case lower({'(Re-)apply (priorly) estimated normalization', 'Apply', 'ReApply'})
            ea_normalize_apply_normalization(options);
        case lower({'FNIRT (Andersson 2010)', 'FNIRT'})
            ea_normalize_fsl(options);
        case lower({'Three-step affine normalization (ANTs; Schonecker 2009)', 'Three-step', 'ThreeStep'})
            ea_normalize_schoenecker(options);
        case lower({'SPM12 DARTEL (Ashburner 2007)', 'SPMDARTEL', 'DARTEL'})
            ea_normalize_spmdartel(options);
        case lower({'SPM12 Segment (Ashburner 2005)', 'SPMSegment', 'Segment'})
            ea_normalize_spmnewseg(options);
        case lower({'SPM12 SHOOT (Ashburner 2011)', 'SPMSHOOT', 'SHOOT'})
            ea_normalize_spmshoot(options);
        otherwise
            warning('Normalization method not recognized...');
            if options.prefs.diary
                diary off;
            end
            return;
    end

    % Compute tone-mapped normalized CT
    if strcmp(options.subj.postopModality, 'CT')
        ea_tonemapct(options, 'norm');
    end

    if options.prefs.diary
        diary off;
    end

    if options.overwriteapproved && isfolder(options.subj.brainshiftDir)
        ea_cprintf('CmdWinWarnings', 'Normalization has been rerun. Please also rerun brain shift correction!\n');
    end
end
