function ea_coregmr(options)
% Wrapper for coregister post-op MRI

% Set fixed anchor image
fixed = options.subj.preopAnat.(options.subj.AnchorModality).coreg;

if ~isfile(fixed)
    warning('No preoperative acquisition found. Coregistration not possible.');
    return
end

% Set moving and output image
moving = cellfun(@(x) options.subj.postopAnat.(x).preproc, fieldnames(options.subj.postopAnat), 'Uni', 0);
out = cellfun(@(x) options.subj.postopAnat.(x).coreg, fieldnames(options.subj.postopAnat), 'Uni', 0);

% Check moving image existence
moving_exists = cellfun(@(x) isfile(x), moving);

% Check registration lock/approval status
out_approved = cellfun(@(x) logical(ea_reglocked(options, x)), out);

% Remove non-existing moving image and approved output image
moving(~moving_exists | out_approved) = [];
out(~moving_exists | out_approved) = [];

% Return if no image remains
if isempty(moving)
    return;
end

doreslice=1;

% Do coregistration
for i = 1:length(moving)
    disp('Coregistering postop MR cor to preop MRI...');
    switch options.coregmr.method
        case 'SPM' % SPM
            % Copy moving image to out image first, since SPM will change
            % the header of the moving image.
            copyfile(moving{i}, out{i});
            ea_coregmr_spm(options, fixed, out{i}, out{i}, doreslice);
        case 'FSL FLIRT' % FSL FLIRT
            ea_coregmr_flirt(options, fixed, moving{i}, out{i});
        case 'FSL BBR' % FSL FLIRT
            ea_coregmr_flirt_bbr(options, fixed, moving{i}, out{i});
        case 'ANTs' % ANTs
            ea_coregmr_ants(options, fixed, moving{i}, out{i}, 0);
        case 'BRAINSFIT' % BRAINSFit
            ea_coregmr_brainsfit(options, fixed, moving{i}, out{i});
        case 'Hybrid SPM & ANTs' % Hybrid SPM -> ANTs
            % Copy moving image to out image first, since SPM will change
            % the header of the moving image.
            copyfile(moving{i}, out{i});
            ea_coregmr_spm(options, fixed, out{i}, out{i}, 0);
            ea_coregmr_ants(options, fixed, out{i}, out{i});
        case 'Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
            % Copy moving image to out image first, since SPM will change
            % the header of the moving image.
            copyfile(moving{i}, out{i});
            ea_coregmr_spm(options, fixed, out{i}, out{i}, 0);
            ea_coregmr_brainsfit(options, fixed, out{i}, out{i});
    end
    disp('Coregistration done.');
end

ea_dumpmethod(options, 'coreg');
