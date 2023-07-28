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

dcm2niix = ea_getExec([basedir, 'dcm2niix'], escapePath = 1);


if strcmp(outdir(end),filesep)
    outdir = outdir(1:end-1);
end

if strcmp(dicomdir(end),filesep)
    dicomdir = dicomdir(1:end-1);
end

cmd=[dcm2niix, ' --progress -f "%f_%p_%z_%t_%s_%d" -i y -b y -v 0 -z y -o ', ea_path_helper(outdir), ' ', ea_path_helper(dicomdir)];

ea_runcmd(cmd);
