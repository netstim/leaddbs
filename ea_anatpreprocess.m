function ea_anatpreprocess(fpth)
% Preprocess of anat image
% Currently includes: reorientation & cropping and bias field correction

%ea_rocrop(fpth); disabled for now.
ea_bias_field_correction(fpth)
