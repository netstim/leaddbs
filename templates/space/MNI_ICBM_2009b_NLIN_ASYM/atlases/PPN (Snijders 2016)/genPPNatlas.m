function genPPNatlas
% About this atlas: Within work of the study Snijders et al. Annals of
% Neurology (2016) entitled "Physiology of Freezing of Gait", the PPN was
% identified by a skilled neuroanatomist based on CHAT positive neuronal
% populations described in Mesulam 1989 ("Human reticular formation: cho-
% linergic neurons of the pedunculopontine and laterodorsal tegmental nuclei
% and some cytochemical comparisons to forebrain cholinergic neurons.").
% This was done directly within ICBM2009b space under additional anatomical
% information provided by the BigBrain 2nd revision template (available
% registered to the 2009c template) and overlay of an HCP template ("T1
% div. by T2" that had been transformed to 2009b space using the MNI T1
% 6thGen NLIN to MNI 2009b NLIN ANTs transform available under
% (https://figshare.com/articles/MNI_T1_6thGen_NLIN_to_MNI_2009b_NLIN_ANTs_transform/3502238).
% Based on this information, the MNI coordinates of x = +/-5.8, y = -33.3,
% z = -20.19 were identified to best represent the PPN within MNI space.
% The resulting location is visualized in the work by Snijders et al., fig.
% 5 panel A. An additional sagittal view of the location is enclosed with
% this atlas. The function below merely creates a sphere of 3 mm radius
% around that MNI location (and since this does not need to be executed
% more than once is written in a very unefficient way).
% - A. Horn 12/2016

nii=ea_load_nii([ea_space,'bb.nii']);
[xx,yy,zz]=ind2sub(size(nii.img),1:numel(nii.img));

XYZ=[xx;yy;zz;ones(1,length(xx))];
XYZ=nii.mat*XYZ;

ix=rangesearch(XYZ(1:3,:)',[5.8, -33.3, -20.19;-5.8, -33.3, -20.19],1);

mkdir('rh')
nii.fname='rh/PPN.nii';
nii.img(:)=0;
nii.img(ix{1})=1;
ea_write_nii(nii);
ea_crop_nii(nii.fname);
mkdir('lh')
nii.fname='lh/PPN.nii';
nii.img(:)=0;
nii.img(ix{2})=1;
ea_write_nii(nii);
ea_crop_nii(nii.fname);
