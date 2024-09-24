function ea_segmask_cleanup(options)
% removes all segmentations if co-registration / normalization were re-done
% this is a "safety" feature
% By Butenko and Li, konstantinmgtu@gmail.com

arguments
    options     % Lead-DBS options for the subject
end

if ispc
    filelist = dir(fullfile(options.subj.subjDir, '**\*.*'));  % get list of files and folders in any subfolder
else
    filelist = dir(fullfile(options.subj.subjDir, '**/*.*'));  % get list of files and folders in any subfolder
end
filelist = filelist(~[filelist.isdir]);  %remove folders from list

% delete all segmask entries from the patient folder
for files_idx = 1:size(filelist,1)
    if contains(filelist(files_idx).name, 'segmask') || contains(filelist(files_idx).name, '-C123_mod') || contains(filelist(files_idx).name, 'label-GM_mask') || contains(filelist(files_idx).name, 'label-WM_mask') || contains(filelist(files_idx).name, 'label-CSF_mask')
        disp("here")
        ea_delete(fullfile(filelist(files_idx).folder,filelist(files_idx).name))
    end
end
