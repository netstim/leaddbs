% extract masked FA maps (mask + 1 -> FA, mask = 0 -> 0)
% save as nifti
% 
% Susanne Schnell
% Medical Physics, Dept. of Radiology, University Hospital Freiburg,
% Germany
%
% Linux
% 25/03/2010
%
% Changes:
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function res = extract_maskedFA

% load tensor file and Roi
[dtdfiles,dir] = uigetfile('*_DTD.mat*','Select all DTDstruct files you want to create the masked FA maps','MultiSelect','On');
if isempty(dtdfiles)
    res = 0;
    return
end

[roifiles,dir2] = uigetfile('*thresh05.mat','Select the masks for masking the FA maps (have to be in same order as the previously given DTDstructs. Only one per DTD file!','MultiSelect','On');

% extract FA and mask it, save as nifti
if iscell(dtdfiles)
    loop = length(dtdfiles);
else
    loop = 1;
end
for m = 1 : loop
    dtd = dtdstruct_read(fullfile(dir,dtdfiles{m}));
    roi = maskstruct_read(fullfile(dir2,roifiles{m}));
    fa = dtdstruct_query(dtd,'getFA');
    fa_AFleft = mrstruct_init('volume', zeros(size(fa.dataAy)),fa);
    inds = find(roi.maskCell{1}(:)==1);
    fa_AFleft.dataAy(inds) = fa.dataAy(inds);
    fa_AFright = mrstruct_init('volume', zeros(size(fa.dataAy)),fa);
    inds = find(roi.maskCell{2}(:)==1);
    fa_AFright.dataAy(inds) = fa.dataAy(inds);
    fa_IOFleft = mrstruct_init('volume', zeros(size(fa.dataAy)),fa);
    inds = find(roi.maskCell{3}(:)==1);
    fa_IOFleft.dataAy(inds) = fa.dataAy(inds);
    fa_IOFright = mrstruct_init('volume', zeros(size(fa.dataAy)),fa);
    inds = find(roi.maskCell{4}(:)==1);
    fa_IOFright.dataAy(inds) = fa.dataAy(inds);
    fa_SLFleft = mrstruct_init('volume', zeros(size(fa.dataAy)),fa);
    inds = find(roi.maskCell{5}(:)==1);
    fa_SLFleft.dataAy(inds) = fa.dataAy(inds);
    fa_SLFright = mrstruct_init('volume', zeros(size(fa.dataAy)),fa);
    inds = find(roi.maskCell{6}(:)==1);
    fa_SLFright.dataAy(inds) = fa.dataAy(inds);
    mrstruct_to_nifti(fa,fullfile(dir,[dtdfiles{m}(1:end-7),'FA.nii']));
    mrstruct_to_nifti(fa_AFleft,fullfile(dir,[dtdfiles{m}(1:end-7),'FA_AFleft.nii']));
    mrstruct_to_nifti(fa_AFright,fullfile(dir,[dtdfiles{m}(1:end-7),'FA_AFright.nii']));
    mrstruct_to_nifti(fa_IOFleft,fullfile(dir,[dtdfiles{m}(1:end-7),'FA_IOFleft.nii']));
    mrstruct_to_nifti(fa_IOFright,fullfile(dir,[dtdfiles{m}(1:end-7),'FA_IOFright.nii']));
    mrstruct_to_nifti(fa_SLFright,fullfile(dir,[dtdfiles{m}(1:end-7),'FA_SLFleft.nii']));
    mrstruct_to_nifti(fa_SLFright,fullfile(dir,[dtdfiles{m}(1:end-7),'FA_SLFright.nii']));
end

res = 1;

