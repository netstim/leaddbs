function ea_rocrop(inputimage)
% Reorient and crop nifti image
% Workaround to fix the 1KB output file issue of dcm2nii

pth = fileparts(inputimage);
if isempty(pth)
    pth = '.';
end

disp(inputimage);

ea_dcm2nii(inputimage, [pth, filesep, 'temp.nii']);

tempf=dir([pth,filesep,'temp.nii']);
tempfsize = tempf.bytes/1024;

% To be safe, retry when the file size smaller than 10 KB
while tempfsize < 10
    disp('dcm2nii returns corrupted file, retrying...')
    ea_dcm2nii(inputimage, [pth, filesep, 'temp.nii']);
    tempf=dir([pth,filesep,'temp.nii']);
    tempfsize = tempf.bytes/1024;
end
movefile([pth,filesep,'temp.nii'],inputimage);
