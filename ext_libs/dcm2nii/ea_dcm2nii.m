function ea_dcm2nii(inputimage)
% Wrapper for dcm2nii, just used for reorientation and cropping currently

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    dcm2nii = [basedir, 'dcm2nii.exe'];
elseif isunix
    dcm2nii = [basedir, 'dcm2nii.', computer];
end

cmd=[dcm2nii, ' -g n -x y ', inputimage];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

[pth,fn,ext]=fileparts(inputimage);
if isempty(pth)
    pth = '.';
end

try
    movefile([pth,filesep,'co',fn,ext],inputimage);
    delete([pth,filesep,'o',fn,ext]); 
catch
    try
        movefile([pth,filesep,'c',fn,ext],inputimage);
    end
end
