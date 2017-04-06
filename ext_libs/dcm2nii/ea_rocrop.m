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

cnt=1; % try ten times, if doesnt work, continue.
failed=0;
% To be safe, retry when the file size smaller than 10 KB
while tempfsize < 10
    disp('dcm2nii returns corrupted file, retrying...')
    ea_dcm2nii(inputimage, [pth, filesep, 'temp.nii']);
    tempf=dir([pth,filesep,'temp.nii']);
    tempfsize = tempf.bytes/1024;
    cnt=cnt+1;
    if cnt>10
        failed=1;
        warning('Giving up. Please check input files, manually.');
        break 
    end
end

if ~failed
    movefile([pth,filesep,'temp.nii'],inputimage);
else
    delete([pth,filesep,'temp.nii']);
end
