function [maskedImage, brainMask] = ea_maskimg(filePath)
% Return path of masked image

% Segment image
filePath = GetFullPath(filePath);
[directory, fileName] = fileparts(filePath);

% Set brain mask and masked image path
if isBIDSFileName(filePath)
    brainMask = setBIDSEntity(filePath, 'mod', parsedStruct.suffix, 'label', 'Brain', 'suffix', 'mask');
    parsedStruct = parseBIDSFilePath(filePath);
    if isfield(parsedStruct, 'acq')
        maskedImage = setBIDSEntity(filePath, 'acq', [], 'label', 'Brain', 'acq', parsedStruct.acq)
    else
        maskedImage = setBIDSEntity(filePath, 'label', 'Brain');
    end
else
    brainMask = [directory, filesep, fileName, '_brainmask.nii'];
    maskedImage = [directory, filesep, fileName, '_brain.nii'];
end

if ~isfile(brainMask)
	ea_genbrainmask(filePath);
end

nii = ea_load_nii(filePath);
mask = ea_load_nii(brainMask);
nii.img = nii.img.*double(mask.img);
nii.fname = maskedImage;
ea_write_nii(nii);
