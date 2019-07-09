function ea_anatpreprocess(fpth)
% Preprocess of anat image.
% Currently includes: reorientation & cropping and bias field correction
%
% USAGE:
%
%    ea_anatpreprocess(fpth)
%
% INPUT:
%    fpth:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ningfei Li, Original file
%       - Daniel Duarte, Documentation

% Cropping and reorientation
ea_rocrop(fpth);

% Bias field correction
ea_bias_field_correction(fpth);
