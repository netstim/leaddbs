function ea_dicm2nii(indir,outdir)
% Wrapper for dicm2nii (Xiangrui Li, https://github.com/xiangruili/dicm2nii)
%
% USAGE:
%
%    ea_dicm2nii(indir,outdir)
%
% INPUTS:
%    indir: DICOM folder
%    outdir: output folder for NIfTI files

dicm2nii(indir, outdir, 0)