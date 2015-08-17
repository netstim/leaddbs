function ea_normalize_fibers(options)
% uses map_coords function by Ged Ridgway (see below)
directory=[options.root,options.patientname,filesep];
[~,preniif]=fileparts(options.prefs.prenii_unnormalized);

if ~exist([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'],'file');
    ea_error('Please run a compatible normalization of the preoperative MRI-volume first. Final (inverse) normalization parameters should be stored as y_ea_inv_normparams.nii inside of the subject folder.');
end

vizz=0; % turn this value to 1 to visualize fiber normalization (option for debugging only, this will drastically slow down the process).

%% check which normalization routine has been used..
% if dartel was used, we need to coregister c2 of b0 and rc2 of anat (since
% deformation fields were estimated for the rc* files and not the native
% anat file.
[options.prefs.b0,options.prefs.prenii_unnormalized]=ea_checkdartelused(options);


% normalize fibers

% get affinematrix from b0 to preop mri
Vfirst=spm_vol([options.root,options.patientname,filesep,options.prefs.b0,',1']);
Vsecond=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized,',1']);
x=spm_coreg(Vfirst,Vsecond);
affinematrix1=Vsecond.mat\spm_matrix(x(:)')*Vfirst.mat;
        %
b0=ea_load_nii([options.root,options.patientname,filesep,options.prefs.b0]);

ysize=size(b0.img,2)+1;

ftr = load([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized]);

reftemplate=[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_1.nii'];
Vmni=spm_vol(reftemplate);
    ysize_mni=Vmni(1).dim(2);

% create (unnormalized) trackvis version
disp('Exporting to TrackVis');
try
reftemplate=[options.root,options.patientname,filesep,options.prefs.b0];
dnii=ea_load_nii(reftemplate);
niisize=size(dnii.img); % get dimensions of reference template.
clear dnii
specs.origin=[0,0,0];
specs.dim=niisize;
try
    H=spm_dicom_headers([root_directory,options.prefs.sampledtidicom]);
    specs.orientation=H{1,1}.ImageOrientationPatient;
catch
    %specs.orientation=[0,1,0,0,0,0];%[0,1,0,-1,0,0];%[1,0,0,0,1,0];
    specs.orientation=[1 0 0 0 -1 0];   %     <----- Original aus example trk_write. Try this one.. %[1,0,0,0,1,0];
end
specs.vox=ftr.vox;
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
%[~,ftrfname]=fileparts(options.prefs.FTR_normalized);
ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
end
disp('Done.');



if vizz
    figure('color','w');
    % plot b0
    subplot(1,3,1);
    title('b0 space');
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img>max(b0.img(:))/7));
    plot3(xx(1:10:end),yy(1:10:end),zz(1:10:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    % plot anat
    anat=ea_load_nii([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
    subplot(1,3,2);
    title('Anat space');
    [xx,yy,zz]=ind2sub(size(anat.img),find(anat.img>max(anat.img(:))/3));
    plot3(xx(1:1000:end),yy(1:1000:end),zz(1:1000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    % plot MNI
    mni=ea_load_nii([options.earoot,'templates',filesep,'mni_hires.nii']);
    subplot(1,3,3);
    title('MNI space');
    [xx,yy,zz]=ind2sub(size(mni.img),find(mni.img>max(mni.img(:))/3));
    % transpose to mm
    XYZ=[xx,yy,zz,ones(length(xx),1)]';
    XYZ=mni.mat*XYZ;
    plot3(XYZ(1,1:10000:end),XYZ(2,1:10000:end),XYZ(3,1:10000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
end

ea_dispercent(0,'Normalizing fibers');
numfibs=length(ftr.curveSegCell);


%ynii=nifti([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);
ynii=nifti([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);
        P = [repmat([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'],3,1),[',1,1';',1,2';',1,3']];
        Vnii = spm_vol(P);
wfibs=cell(length(ftr.curveSegCell),1);
for fib=1:numfibs
    
    ea_dispercent(fib/numfibs);
    
    %% transpose from freiburg to spm notation.
    wfibs{fib}=[ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,2),ftr.curveSegCell{fib}(:,3),ones(length(ftr.curveSegCell{fib}),1)];
    %wfibs{fib}=[ftr.curveSegCell{fib}(:,2),ysize-ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,3),ones(length(ftr.curveSegCell{fib}),1)];
    if vizz
       thisfib=wfibs{fib}';
        subplot(1,3,1)
        plot3(thisfib(1,:),thisfib(2,:),thisfib(3,:),'-','color',[0.1707    0.2919    0.7792]);
    end
    
    %% first apply affine transform from b0 to prenii
    wfibs{fib}=affinematrix1*wfibs{fib}';
    if vizz
        thisfib=wfibs{fib};
        subplot(1,3,2)
        plot3(thisfib(1,:),thisfib(2,:),thisfib(3,:),'-','color',[0.1707    0.2919    0.7792]);
    end
    %% -> coordinates are now in voxel-space of single subject anat file.
    
    
 
 
    %% map from prenii voxelspace to mni millimeter space   
    
    wfibs{fib} = vox2mm_mni(wfibs{fib},Vnii,ynii)';
    
    
    
    if vizz
        thisfib=wfibs{fib}';
        subplot(1,3,3)
        plot3(thisfib(1,:),thisfib(2,:),thisfib(3,:),'-','color',[0.1707    0.2919    0.7792]);
    end
    
    
    %% map from mni millimeter space to mni voxel space (only needed for trackvis convertion).
    wfibsvox{fib}=[wfibs{fib},ones(size(wfibs{fib},1),1)]';
    wfibsvox{fib}=Vmni(1).mat\wfibsvox{fib};
    wfibsvox{fib}=wfibsvox{fib}(1:3,:)';
    wfibsvox{fib}=[wfibsvox{fib}(:,1),ysize_mni-wfibsvox{fib}(:,2),wfibsvox{fib}(:,3)];
    
    
    %% cleanup
    wfibs{fib}=wfibs{fib}(:,1:3);
   if vizz; drawnow; end
end
ea_dispercent(100,'end');


wfibs=wfibs';
wfibsvox=wfibsvox';
normalized_fibers_mm=wfibs; clear wfibs
normalized_fibers_vox=wfibsvox; clear wfibsvox

save([options.root,options.patientname,filesep,options.prefs.FTR_normalized],'normalized_fibers_mm');
save([options.root,options.patientname,filesep,'vox_',options.prefs.FTR_normalized],'normalized_fibers_vox');


 

% create trackvis version
try
reftemplate=[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_1.nii'];
dnii=ea_load_nii(reftemplate);
niisize=size(dnii(1).img); % get dimensions of reference template.
clear dnii
specs.origin=[0,0,0];
specs.dim=niisize;
try
    H=spm_dicom_headers([root_directory,options.prefs.sampledtidicom]);
    specs.orientation=H{1,1}.ImageOrientationPatient;
catch
    specs.orientation=[0,1,0,0,0,0]; %[0,1,0,-1,0,0];%;[0,1,0,-1,0,0] [0,1,0,0,0,0];
    specs.orientation=[1,0,0,0,1,0];
    specs.orientation=[1,0,0,1,0,0];
    specs.orientation=[1,0,0,0,1,0];
    %specs.orientation=[1 0 0 0 -1 0];   %     <----- Original aus example trk_write. Try this one.. %[1,0,0,0,1,0];
end
specs.vox=ftr.vox;

[~,ftrfname]=fileparts(options.prefs.FTR_normalized);
ea_ftr2trk(['vox_',ftrfname],directory,specs,options); % export normalized ftr to .trk

movefile([directory,'vox_',ftrfname,'.trk'],[directory,ftrfname,'.trk']);

end
delete([options.root,options.patientname,filesep,'vox_',options.prefs.FTR_normalized]);

disp('Done.');





function [useb0,useanat]=ea_checkdartelused(options)
directory=[options.root,options.patientname,filesep];

dartelused=0;
try
    load([directory,'ea_normmethod_applied']);
    
    if strcmp(norm_method_applied{end},'ea_normalize_spmdartel')
        dartelused=1;
    end
end
if dartelused
    % segment b0.
    if ~exist([directory,'c2',options.prefs.b0],'file');
        disp('Segmenting B0 file for DARTEL import space coregistration...');
        ea_newseg(directory,options.prefs.b0,0,options);
        delete([directory,'c4',options.prefs.b0]);
        delete([directory,'c5',options.prefs.b0]);
        disp('Done.');
    end
    
end

if ~dartelused
    useb0=options.prefs.b0;
    useanat=options.prefs.prenii_unnormalized;
else
    useb0=['c2',options.prefs.b0];
    useanat=['rc2',options.prefs.prenii_unnormalized];
end


function coord = vdox2mm_mni(coord,Vnii,ynii)



ixs = double(coord(1:3, :));

% old method
for i = 1:3
coord(i,:)=ynii.dat(sub2ind(size(ynii.dat),ixs(1,:)',ixs(2,:)',ixs(3,:)',ones(size(ixs,2),1),repmat(i,size(ixs,2),1)));
end

function coord = vox2mm_mni(coord, Vnii,ynii)
% new method
ixs_new = double(coord(1:3, :));
coord=[spm_sample_vol(Vnii(1),ixs_new(1,:),ixs_new(2,:),ixs_new(3,:),1);
    spm_sample_vol(Vnii(2),ixs_new(1,:),ixs_new(2,:),ixs_new(3,:),1);
    spm_sample_vol(Vnii(3),ixs_new(1,:),ixs_new(2,:),ixs_new(3,:),1)];


% %old method
% ixs_old = round(coord(1:3, :));
% coord=zeros(3,size(coord,2));
% for i = 1:3
% coord(i,:)=ynii.dat(sub2ind(size(ynii.dat),ixs_old(1,:)',ixs_old(2,:)',ixs_old(3,:)',ones(size(ixs_old,2),1),repmat(i,size(ixs_old,2),1)));
% end
  
