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
    dcm2niix = ea_path_helper([basedir, 'dcm2niix.exe']);
else
    dcm2niix = [basedir, 'dcm2niix.', computer('arch')];
end

if strcmp(outdir(end),filesep)
    outdir = outdir(1:end-1);
end
if strcmp(dicomdir(end),filesep)
    dicomdir = dicomdir(1:end-1);
end
cmd=[dcm2niix, ' -f "%f_%p_%z_%t_%s_%d" -i y -b y -v 0 -z y y', ' -o ', ea_path_helper(outdir), ' ', ea_path_helper(dicomdir)];

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
