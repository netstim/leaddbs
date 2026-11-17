function [faNative, faCoreg] = ea_getfapaths(options)
% Return absolute paths of the FA map in DWI space (faNative) and of the FA
% map coregistered to the anchor (faCoreg). Both are '' when not defined.
%
% The coregistered FA must always share the B0->anchor transform, so every
% function that reads or writes it should use these paths.

faNative = '';
faCoreg = '';

if ~isfield(options, 'subj') || ~isfield(options.subj, 'subjDir')
    return;
end
subjDir = options.subj.subjDir;

if isfield(options, 'prefs') && isfield(options.prefs, 'fa') && ~isempty(options.prefs.fa)
    faNative = makeAbsolute(options.prefs.fa, subjDir);
end

if isfield(options, 'prefs') && isfield(options.prefs, 'fa2anat') && ...
        ~isempty(options.prefs.fa2anat) && ~strcmp(options.prefs.fa2anat, 'fa2anat.nii')
    faCoreg = makeAbsolute(options.prefs.fa2anat, subjDir);
elseif isfield(options.subj, 'subjId')
    % BIDS default (same name as written by ea_ensure_fa_and_fa2anat)
    faCoreg = fullfile(subjDir, 'coregistration', 'anat', ...
        ['sub-', options.subj.subjId, '_ses-preop_space-anchorNative_fa.nii']);
end


function p = makeAbsolute(p, subjDir)
isAbsolute = startsWith(p, filesep) || ...
    ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once')) || startsWith(p, '\\');
if ~isAbsolute
    p = fullfile(subjDir, p);
end
