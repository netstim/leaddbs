% give the size of several rois in one folder

% Susanne Schnell
% Medical Physics
% Dept. of Radiology, University Hospital Freiburg, Germany

% 18/03/2010
% Linux
% Changes: 

function size_roi = roi_size(thresh)

if nargin == 0
    thresh = 0;
end

[files,dir] = uigetfile({'*.nii';'*.img'},'Select the rois to evaluate the size of it!','MultiSelect','on');
if isempty(files)
    return
end
if iscell(files)
    numfor = length(files);   
else
    numfor = 1;
end
for m = 1 : numfor
    roi = nifti_to_mrstruct('volume',fullfile(dir, files{m}));
    size_roi.name{m,:} = files{m};
    inds = find(roi.dataAy(:) > thresh);
    roi.dataAy(:) = zeros(size(roi.dataAy(:)));
    roi.dataAy(inds) = 1;
    size_roi.size(m) = sum(roi.dataAy(inds));
end

