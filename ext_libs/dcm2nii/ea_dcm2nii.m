function ea_dcm2nii(inputimage)
% Wrapper for dcm2nii, just used for reorientation and cropping currently

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    dcm2nii = [basedir, 'dcm2nii.exe'];
else
    dcm2nii = [basedir, 'dcm2nii.', computer('arch')];
end

cmd=[dcm2nii, ' -g n -x y ', inputimage];

fprintf('\nReorient and crop image...\n')
if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

[pth,fn,ext]=fileparts(inputimage);
if isempty(pth)
    pth = '.';
end

if exist([pth,filesep,'co',fn,ext], 'file')
    movefile([pth,filesep,'co',fn,ext],inputimage);
    delete([pth,filesep,'o',fn,ext]);
elseif exist([pth,filesep,'c',fn,ext], 'file')
    movefile([pth,filesep,'c',fn,ext],inputimage);
else
	disp('No need to cropping!');
end
