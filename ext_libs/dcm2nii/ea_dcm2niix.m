function ea_dcm2niix(inputimage,ofolder)
% Wrapper for dcm2niix

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    dcm2nii = [basedir, 'dcm2niix.exe'];
elseif isunix
    dcm2nii = [basedir, 'dcm2niix.', computer];
end

cmd=[dcm2nii, ' -o ',ofolder,' -z n -x y ',inputimage];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end


% -o option seems not to work in dcm2niix yet. Move files manually:
di=dir([inputimage,'*.nii']);
for d=1:length(di)
    movefile([inputimage,di(d).name],[ofolder,di(d).name]);
end

di=dir([inputimage,'*.bval']);
for d=1:length(di)
    movefile([inputimage,di(d).name],[ofolder,di(d).name]);
end

di=dir([inputimage,'*.bvec']);
for d=1:length(di)
    movefile([inputimage,di(d).name],[ofolder,di(d).name]);
end

% try
%     delete([pth,filesep,'o',fn,ext]);
%     movefile([pth,filesep,'co',fn,ext],inputimage);
% catch
%     try
%         movefile([pth,filesep,'c',fn,ext],inputimage);
%     end
% end
