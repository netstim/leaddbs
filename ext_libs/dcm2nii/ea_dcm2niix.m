function ea_dcm2niix(dicomdir, outdir)
% Wrapper for dcm2niix

if nargin < 2
    outdir = dicomdir;
end

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

ea_libs_helper;

basedir = [fileparts(mfilename('fullpath')), filesep];

if ispc
    dcm2niix = [basedir, 'dcm2niix.exe'];
else
    dcm2niix = [basedir, 'dcm2niix.', computer('arch')];
end

cmd=[dcm2niix, ' -o ', ea_path_helper(outdir), ' -z n -x y ', ea_path_helper(dicomdir)];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
