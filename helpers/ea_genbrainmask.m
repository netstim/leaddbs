function ea_genbrainmask(filePath)
% Generate brain mask based on SPM New Segment

% Segment image
filePath = GetFullPath(filePath);
[directory, fileName] = fileparts(filePath);
ea_newseg(filePath, 0, 1);

% Load segmentations
c1 = ea_load_nii([directory, filesep, 'c1', fileName, '.nii']);
c2 = ea_load_nii([directory, filesep, 'c2', fileName, '.nii']);
c3 = ea_load_nii([directory, filesep, 'c3', fileName, '.nii']);

% Join masks
c1.img = c1.img + c2.img + c3.img;
c1.img = c1.img > 0.6;

% Set header accordingly
c1.dt = [2 0]; % unit8 according to spm_type
c1.descrip = 'Brain Mask';
c1.pinfo(1:2) = [1,0]; % uint8 is enough for output values, no need for scaling

% Set output file name
if isBIDSFileName(filePath)
    parsedStruct = parseBIDSFilePath(filePath);
    c1.fname = setBIDSEntity(filePath, 'mod', parsedStruct.suffix, 'label', 'Brain', 'suffix', 'mask');
    movefile([directory, filesep, 'c1', fileName, '.nii'], setBIDSEntity(filePath, 'mod', parsedStruct.suffix, 'label', 'GM', 'suffix', 'mask'));
    movefile([directory, filesep, 'c2', fileName, '.nii'], setBIDSEntity(filePath, 'mod', parsedStruct.suffix, 'label', 'WM', 'suffix', 'mask'));
    movefile([directory, filesep, 'c3', fileName, '.nii'], setBIDSEntity(filePath, 'mod', parsedStruct.suffix, 'label', 'CSF', 'suffix', 'mask'));
else
    c1.fname = [directory, filesep, fileName, '_brainmask.nii'];
end

% Write image
ea_write_nii(c1);
