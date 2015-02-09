function varargout=ea_ft_deterministictracking_kroon(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Deterministic Fibertracking (Kroon).';
    varargout{2}={'SPM8','SPM12'};
    return
end

%% mask for tracking


%% mask for tracking

directory=[options.root,options.patientname,filesep];


% 'new segment' options.prefs.prenii_unnormalized
matlabbatch{1}.spm.tools.preproc8.channel.vols = {[directory,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spm('dir'),filesep,'toolbox/Seg/TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 1];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spm('dir'),filesep,'toolbox/Seg/TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 1];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spm('dir'),filesep,'toolbox/Seg/TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spm('dir'),filesep,'toolbox/Seg/TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spm('dir'),filesep,'toolbox/Seg/TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spm('dir'),filesep,'toolbox/Seg/TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];


jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;


%% Coreg options.prefs.prenii_unnormalized to b0 (for label.mat and FTR-Normalization)


matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[directory,'c1',options.prefs.prenii_unnormalized,',1']
    [directory,'c2',options.prefs.prenii_unnormalized,',1']
    [directory,'c3',options.prefs.prenii_unnormalized,',1']
    [directory,'c4',options.prefs.prenii_unnormalized,',1']
    [directory,'c5',options.prefs.prenii_unnormalized,',1']
    };
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'rb0';

jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs;



%% prepare DTI
Vdti=spm_vol([options.root,options.patientname,filesep,options.prefs.dti]);
Xdti=spm_read_vols(Vdti);

bval=load([options.root,options.patientname,filesep,options.prefs.bval]);
bvec=load([options.root,options.patientname,filesep,options.prefs.bvec]);

for i=1:size(Xdti,4)
   DTIdata(i).VoxelData=single(squeeze(Xdti(:,:,:,i)));
   DTIdata(i).Gradient=bvec(:,i);
   DTIdata(i).Bvalue=bval(:,i);
end


% Constants DTI
parametersDTI=[];
parametersDTI.BackgroundTreshold=150;
parametersDTI.WhiteMatterExtractionThreshold=0.10;
parametersDTI.textdisplay=true;

% Perform DTI calculation
[ADC,FA,VectorF,DifT]=DTI(DTIdata,parametersDTI);


% load mask
Vmask=spm_vol([directory,'rb0c2',options.prefs.prenii_unnormalized,',1']);
Xmask=spm_read_vols(Vmask);
Xmask=Xmask>0.5;
Xmask(1,:,:)=0;
Xmask(end,:,:)=0;
Xmask(:,1,:)=0;
Xmask(:,end,:)=0;
Xmask(:,:,1)=0;
Xmask(:,:,end)=0;

% Fiber Tracking Constants
parametersFT=[];
parametersFT.FiberLengthMax=600;
parametersFT.FiberLengthMin=6;
parametersFT.DeviationAngleMax=1;
parametersFT.Step=0.4;
parametersFT.FiberTrackingComputationTreshold=0.125;
parametersFT.Sampling=2;
parametersFT.textdisplay=true;
parametersFT.FAmin=0.2;
parametersFT.FAmax=0.8;

% Perform FT
fibers=ea_FT(FA,VectorF,Xmask,parametersFT);


save([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized],'fibers');




%% export .trk copy for trackvis visualization

dnii=load_nii([directory,options.prefs.b0]);
niisize=size(dnii.img); % get dimensions of reference template.
specs.origin=[0,0,0];
specs.dim=niisize;
specs.vox=dnii.hdr.dime.pixdim(2:4);
try
    H=spm_dicom_headers([directory,prefs.sampledtidicom]);
    specs.orientation=H{1,1}.ImageOrientationPatient;
catch
    specs.orientation=[1,0,0,0,1,0];
end
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
disp('Done.');





