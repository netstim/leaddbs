function [coords,trajvector,trajectory,tramat]=ea_reconstruct(patientname,options,side)
% This function is the heart of Lead-DBS, reconstructing the electrode
% lead trajectory for one side (left/right) of the MR-data. It reads two MR
% images from a folder called 'patientname' and iteratively reconstructs a
% line in 3D-space that best describes the electrode trajectory.
% __________________________________________________________________________________
%
% Inputs:   patientname     – String of folder and root of filenames. Files
%                             should be called e.g.
%                             'MustermannMax/MustermannMax_tra_brain_A3_final.nii'
%                             and
%                             'MustermannMax/MustermannMax_cor_brain_A3_final_opt.nii'
%                             and reside within the folder specified by
%                             options.root. Note that an exact
%                             normalization into MNI-space is crucial for
%                             Lead to work correctly.
%           options         – Struct containing various options, see e.g.
%                             ea_defaultoptions.m
%           side            – which side of the brain shall be
%                             reconstructed. 1 > right hemisphere, 2 > left hemisphere. 
% ----------------------------------------------------------------------------------
% 
% Outputs:  coords          – 8x3 vector of electrode coordinates in
%                             mm-representations within MNI-space (if
%                             MR-images have been normalized correctly).
%           trajvector      – 3 element vector describing the traversing
%                             direction of the lead trajectory.
%           trajectory         – nx3 vector describing the fitted line of the
%                             trajectory.
%           tramat          – 4x4 matrix describing the normalization of
%                             MR-images. This is used to reconstruct the
%                             distances of the electrode contacts in
%                             MNI-space (since the distance has changed
%                             from e.g. 2mm in native space due to
%                             normalization).
%                             
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


tra_nii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
try
cor_nii=load_untouch_nii([options.root,options.prefs.patientdir,filesep,options.prefs.cornii]);
end

imat=zeros([size(tra_nii.img,1),size(tra_nii.img,2),size(tra_nii.img,3),2]);
imat(:,:,:,1)=tra_nii.img;
try
imat(:,:,:,2)=cor_nii.img;
end

ea_showdis('Preparing contrasted volume...',options.verbose);

tra_nii.img=ea_gencontrastimage(imat,options.axiscontrast);

trajectory=[]; % empty initialization.
for refine=0:options.refinesteps
[trajectory,trajvector]=ea_reconstruct_trajectory(trajectory,tra_nii,side,refine,options);
end






%% determine height of last electrode

%detdiams=detrend(diams);


if options.verbose>1; di=figure('name','Finding local maxima in diameters...','numbertitle','off'); end
if options.verbose>2; close(di); end





% find local maxima in diameters.


% first, calculate distance between contacts.
% rename matfile to text
try
tramat=load([options.root,patientname,filesep,'ea_normparams']);
tramat=tramat.M;
catch
    try
    movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.mat'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt']);
% load matrices
tramat=load([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt']);
% rename the file to .mat again
movefile([options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.txt'],[options.root,patientname,filesep,patientname,'_tra_brain_A3_mo_final.mat']);
    catch
    try
        
       tramat=load([options.root,patientname,filesep,patientname,'mmflirt_transform.txt']);
 
    catch
    tramat=eye(4);
    end
    end
end



[~,~,dist]=ea_calc_distance(options.elspec.eldist,trajvector,tramat(1:3,1:3),[options.root,patientname,filesep,options.prefs.tranii]);
% zdist is the distance between electrodes in z-direction.
zdist=dist/norm(trajvector);




% 
% % transform trajectory to mm space:
% 
%             if ~isempty(trajectory)
%                 trajectory=ea_map_coords(trajectory', [options.root,patientname,filesep,options.prefs.tranii])';
%             end

[coords,goodz]=ea_sample_cuboid(trajectory,trajvector,options);




% determine coords by goodz, trajectory and correction term.
% 
% correction=[0,0,0];
% try
%     coords=ea_findcoords(goodz(1),trajectory,trajvector,dist,correction,options);
% catch
%     ea_showdis('Coords not found.',options.verbose);
%     return
% end





