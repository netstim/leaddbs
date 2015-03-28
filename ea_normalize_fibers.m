function ea_normalize_fibers(options)
% uses map_coords function by Ged Ridgway (see below)

[~,preniif]=fileparts(options.prefs.prenii_unnormalized);
if ~exist([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'],'file');
    ea_error('Please run a compatible normalization of the preoperative MRI-volume first. Final (inverse) normalization parameters should be stored as y_ea_inv_normparams.nii inside of the subject folder.');
end



% normalize fibers

% get affinematrix from b0 to preop mri
Vfirst=spm_vol([options.root,options.patientname,filesep,options.prefs.b0,',1']);
Vsecond=spm_vol([options.root,options.patientname,filesep,options.prefs.prenii_unnormalized,',1']);
x=spm_coreg(Vfirst,Vsecond);
affinematrix1=Vsecond.mat\spm_matrix(x(:)')*Vfirst.mat;
        %


b0=load_untouch_nii([options.root,options.patientname,filesep,options.prefs.b0]);

ysize=size(b0.img,2)+1;

ftr = load([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized]);
ea_dispercent(0,'Normalizing fibers');
numfibs=length(ftr.curveSegCell);




for fib=1:numfibs
    
    ea_dispercent(fib/numfibs);
    
    %% transpose from freiburg to spm notation.
    
    wfibs{fib}=[ftr.curveSegCell{fib}(:,2),ysize-ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,3),ones(length(ftr.curveSegCell{fib}),1)];
 
    
    
    %% first apply affine transform from b0 to prenii
    wfibs{fib}=[affinematrix1*wfibs{fib}'];
    
    %% -> coordinates are now in voxel-space of single subject prenii file.
    
 
 
    %% map from prenii voxelspace to mni-millimeter space   
    wfibs{fib} = ea_map_coords(wfibs{fib}, '', [options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'])';
    
    %% cleanup
    wfibs{fib}=wfibs{fib}(:,1:3);
    
end
ea_dispercent(100,'end');


wfibs=wfibs';
normalized_fibers_mm=wfibs; clear wfibs
save([options.root,options.patientname,filesep,options.prefs.FTR_normalized],'normalized_fibers_mm');


 

% create trackvis version
reftemplate=[spm('dir'),filesep,'canonical',filesep,'single_subj_T1.nii'];
dnii=load_nii(reftemplate);
niisize=size(dnii.img); % get dimensions of reference template.
clear dnii
specs.origin=[0,0,0];
specs.dim=niisize;
try
    H=spm_dicom_headers([root_directory,options.prefs.sampledtidicom]);
    specs.orientation=H{1,1}.ImageOrientationPatient;
catch
    specs.orientation=[1,0,0,0,1,0];
end
[~,ftrfname]=fileparts(options.prefs.FTR_normalized);
ea_ftr2trk(ftrfname,directory,specs); % export normalized ftr to .trk

ftr2trk(['w',ftrfilename,usedartel],directory,['w',ftrfilename,usedartel],specs); % export normalized ftr to .trk
disp('Done.');

