function ea_segmask_cleanup(options)
% removes all segmentations if co-registration / normalization were re-done
% this is a "safety" feature
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for the subject
end

pattern = '(segmask|-C123_mod|label-(GM|WM|CSF)_mask)';
maskFiles = ea_regexpdir(options.subj.subjDir, pattern, 1, 'f');
if ~isempty(maskFiles)
    ea_delete(maskFiles);
    ea_cprintf('Comments', 'Segmask cleaned up for "%s".\n', options.subj.subjId);
end