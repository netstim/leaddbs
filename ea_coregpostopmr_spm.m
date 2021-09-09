function varargout = ea_coregpostopmr_spm(options, fixed, moving, out, doreslice)
% Wrapper for SPM registration of post-op MRI

if ischar(options) % return name of method.
    varargout{1} = 'SPM (Friston 2007)';
    return
end

% Available cost functions
% 'mi'  - Mutual Information
% 'nmi' - Normalised Mutual Information
% 'ecc' - Entropy Correlation Coefficient
% 'ncc' - Normalised Cross Correlation
costfuns = {'mi', 'nmi'};

for i=1:length(costfuns)-1
    % Estimate without reslice
    ea_spm_coreg(options, moving, fixed, ...
        costfuns{i}, 0, {''}, options.prefs.mrcoreg.writeoutcoreg);
    disp(['*** Coregistration (Estimate) completed with cost function ''',costfuns{i},'''.']);
end

% Estimate and reslice if necessary
ea_spm_coreg(options, moving, fixed, ...
    costfuns{end}, doreslice, {''}, options.prefs.mrcoreg.writeoutcoreg);

% Verbose and move file
if doreslice
    spmOutput = regexprep(moving, ['([^\', filesep, ']+\.nii(\.gz)?$)'], 'r$1');
    if isfile(spmOutput) && ~strcmp(spmOutput, out)
        movefile(spmOutput, out);
    end
    disp(['*** Coregistration (Estimate and Reslice) completed with cost function ''',costfuns{end},'''.']);
else
    disp(['*** Coregistration (Estimate) completed with cost function ''', costfuns{end},'''.']);
end
