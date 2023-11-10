function ea_convert_atlas2segmask(gm_mask, segmask_file, threshold)

% convert probabilistic GM mask to segmaks for OSS-DBS
% IMPORTANT: all other tissue is assumed WM, no CSF!

prob_gm_mask = ea_load_nii(gm_mask);

% initialize
Segment_GM_WM_CSF = prob_gm_mask;

% change the name
Segment_GM_WM_CSF.fname = segmask_file;

data = Segment_GM_WM_CSF.img;

% default WM except background
data(Segment_GM_WM_CSF.img(:) > 0) = 2;

% assign GM
data(Segment_GM_WM_CSF.img(:) >= threshold) = 1;

Segment_GM_WM_CSF.img = data;

%save
ea_write_nii(Segment_GM_WM_CSF);