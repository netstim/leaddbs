function ea_normalize(options)
% Entry function to run normalization

if ~ea_reglocked(options, options.subj.preopAnat.(options.subj.AnchorModality).norm)
    % Setup log
    ea_mkdir(fileparts(options.subj.norm.log.logBaseName));
    diary([options.subj.norm.log.logBaseName, datestr(now, 'yyyymmddTHHMMss'), '.log']);
    
    % Do coregistration
    switch lower(options.normalize.method)
        case lower({'ANTs (Avants 2008)', 'ANTs'})
            ea_normalize_ants(options);
        case lower({'(Re-)apply (priorly) estimated normalization.', 'Apply', 'ReApply'})
            ea_normalize_apply_normalization(options);
        case lower({'FNIRT (Andersson 2010)', 'FNIRT'})
            ea_normalize_fsl(options);
        case lower({'MAGeT Brain-like Normalization (Chakravarty 2013)', 'MAGeTNorm'})
            ea_normalize_maget(options);
        case lower({'MAGeT Brain-like Segmentation/Normalization DISTAL atlas (Chakravarty 2013, Ewert 2016)', 'MAGeTSeg'})
            ea_normalize_maget_segment(options);
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
            diary off;
            return;
    end

    % Dump method
    ea_dumpmethod(options, 'norm');

    % Compute tone-mapped normalized CT
    if options.modality == 2
        ea_tonemapct(options, 'norm');
    end

    diary off;
end
