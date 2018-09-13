function ea_anatpreprocess(fpth)
% Preprocess of anat image
% Currently includes: reorientation & cropping and bias field correction

% Cropping and reorientation
ea_rocrop(fpth);

% Bias field correction
ea_bias_field_correction(fpth);
