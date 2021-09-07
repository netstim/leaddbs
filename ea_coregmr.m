function ea_coregmr(options)
% wrapper for coreg routines

% in CT imaging, coregistration is done elsewhere.
% also ignore when there is no tra/cor/sag existing (normal conn study)
if options.modality == 2 || ~isfield(options.prefs,'tranii_unnormalized')
    return
end

fixed = options.subj.preopAnat.(options.subj.AnchorModality).coreg;

if ~isfile(fixed)
    warning('No preoperative acquisition found. Coregistration not possible.');
    return
end
    
moving = cellfun(@(x) options.subj.postopAnat.(x).preproc, fieldnames(options.subj.postopAnat), 'UniformOutput', false);
out = cellfun(@(x) options.subj.postopAnat.(x).coreg, fieldnames(options.subj.postopAnat), 'UniformOutput', false);

moving_exists = cellfun(@(x) isfile(x), moving);
out_approved = cellfun(@(x) logical(ea_coreglocked(options, x)), out);

moving(~moving_exists | out_approved) = [];
out(~moving_exists | out_approved) = [];

if isempty(moving); return; end

doreslice=1;

for i = 1:length(moving)
    disp('Coregistering postop MR cor to preop MRI...');
    switch options.coregmr.method
        case 'SPM' % SPM
            ea_coregmr_spm(options, fixed, moving{i}, out{i}, doreslice);
        case 'FSL FLIRT' % FSL FLIRT
            ea_coregmr_flirt(options, fixed, moving{i}, out{i});
        case 'FSL BBR' % FSL FLIRT
            ea_coregmr_flirt_bbr(options, fixed, moving{i}, out{i});
        case 'ANTs' % ANTs
            ea_coregmr_ants(options, fixed, moving{i}, out{i}, 0);
        case 'BRAINSFIT' % BRAINSFit
            ea_coregmr_brainsfit(options, fixed, moving{i}, out{i});
        case 'Hybrid SPM & ANTs' % Hybrid SPM -> ANTs
            ea_coregmr_spm(options, fixed, moving{i}, out{i}, 0); % dont use doreslice here to refrain for doing two interpolations.
            ea_coregmr_ants(options, fixed, moving{i}, out{i});
        case 'Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
            ea_coregmr_spm(options, fixed, moving{i}, out{i}, 0); % dont use doreslice here to refrain for doing two interpolations.
            ea_coregmr_brainsfit(options, fixed, moving{i}, out{i});
    end
    disp('Coregistration done.');
end

ea_dumpmethod(options, 'coreg');
