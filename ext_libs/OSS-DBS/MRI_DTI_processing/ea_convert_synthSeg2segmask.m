function ea_convert_synthSeg2segmask(SynthSeg_segmask_image, segmask_output)

% convert SynthSeg segmentations to GM/WM/CSF
% IMPORTANT: lesions are not segmented, check for atrophies beforehand

SynthSeg_segmask = ea_load_nii(SynthSeg_segmask_image);

% initialize
Segment_GM_WM_CSF = SynthSeg_segmask;

% change the name
Segment_GM_WM_CSF.fname = segmask_output;

data = Segment_GM_WM_CSF.img;

% see SynthSeg labels at https://surfer.nmr.mgh.harvard.edu/fswiki/SynthSeg

% default GM except background
data(Segment_GM_WM_CSF.img(:) > 0) = 1;

% WM
data(Segment_GM_WM_CSF.img(:) == 2) = 2;  % left cerebral WM
data(Segment_GM_WM_CSF.img(:) == 7) = 2;  % left cerebellar WM
data(Segment_GM_WM_CSF.img(:) == 41) = 2; % right cerebral WM
data(Segment_GM_WM_CSF.img(:) == 46) = 2; % right cerebellar WM
data(Segment_GM_WM_CSF.img(:) == 16) = 2; % brainstem

% assign ventricles and CSF to CSF
data(Segment_GM_WM_CSF.img(:) == 4) = 3;
data(Segment_GM_WM_CSF.img(:) == 5) = 3;
data(Segment_GM_WM_CSF.img(:) == 14) = 3;
data(Segment_GM_WM_CSF.img(:) == 15) = 3;
data(Segment_GM_WM_CSF.img(:) == 24) = 3;
data(Segment_GM_WM_CSF.img(:) == 43) = 3;
data(Segment_GM_WM_CSF.img(:) == 44) = 3;

Segment_GM_WM_CSF.img = data;

ea_write_nii(Segment_GM_WM_CSF);