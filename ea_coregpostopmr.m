function ea_coregpostopmr(options)
% Entry function to coregister post-op MRI to pre-op MRI

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

% Setup log
ea_mkdir(fileparts(options.subj.coreg.log.logBaseName));
diary([options.subj.coreg.log.logBaseName, 'MR', datestr(now, 'yyyymmddTHHMMss'), '.log']);

% Do coregistration
for i = 1:length(moving)
    disp('Coregistering postop MR cor to preop MRI...');
    switch options.coregmr.method
        case 'SPM' % SPM
            % Copy moving image to out image first, since SPM will change
            % the header of the moving image.
            copyfile(moving{i}, out{i});
            ea_coregpostopmr_spm(options, fixed, out{i}, out{i}, doreslice);
        case 'FSL FLIRT' % FSL FLIRT
            ea_coregpostopmr_flirt(options, fixed, moving{i}, out{i});
        case 'FSL BBR' % FSL FLIRT
            ea_coregpostopmr_flirtbbr(options, fixed, moving{i}, out{i});
        case 'ANTs' % ANTs
            ea_coregpostopmr_ants(options, fixed, moving{i}, out{i}, 0);
        case 'BRAINSFIT' % BRAINSFit
            ea_coregpostopmr_brainsfit(options, fixed, moving{i}, out{i});
        case 'Hybrid SPM & ANTs' % Hybrid SPM -> ANTs
            % Copy moving image to out image first, since SPM will change
            % the header of the moving image.
            copyfile(moving{i}, out{i});
            ea_coregpostopmr_spm(options, fixed, out{i}, out{i}, 0);
            ea_coregpostopmr_ants(options, fixed, out{i}, out{i});
        case 'Hybrid SPM & BRAINSFIT' % Hybrid SPM -> Brainsfit
            % Copy moving image to out image first, since SPM will change
            % the header of the moving image.
            copyfile(moving{i}, out{i});
            ea_coregpostopmr_spm(options, fixed, out{i}, out{i}, 0);
            ea_coregpostopmr_brainsfit(options, fixed, out{i}, out{i});
    end
    disp('Coregistration done.');
end

ea_dumpmethod(options, 'coreg');

diary off;
