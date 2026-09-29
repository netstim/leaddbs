function archived = ea_cleanup_normtransforms(options)
% Retire inactive normalization fields only after the active pair exists.
active = ea_gettransformfiles(options);
if ~isfile(active.forward) || ~isfile(active.inverse)
    error('LeadDBS:MissingNormalization', 'Cannot retire old warps before the current pair exists.');
end
bases = {options.subj.norm.transform.forwardBaseName, options.subj.norm.transform.inverseBaseName};
suffixes = {'ants.nii.gz','ants.mat','spm.nii','fnirt.nii.gz'};
candidates = {};
for b = 1:numel(bases)
    for s = 1:numel(suffixes)
        candidates{end+1} = [bases{b}, suffixes{s}]; %#ok<AGROW>
    end
end
folder = fileparts(active.forward);
candidates = [candidates, {fullfile(folder,'y_ea_normparams.nii'), ...
    fullfile(folder,'y_ea_inv_normparams.nii')}];
archived = {};
archive = '';
for i = 1:numel(candidates)
    old = candidates{i};
    if isfile(old) && ~strcmp(old,active.forward) && ~strcmp(old,active.inverse)
        if isempty(archive)
            logFolder = fileparts(options.subj.norm.log.method);
            archive = tempname(logFolder);
            mkdir(archive);
        end
        [~,name,ext] = fileparts(old);
        destination = fullfile(archive,[name,ext]);
        movefile(old,destination);
        archived{end+1} = destination; %#ok<AGROW>
    end
end
end

