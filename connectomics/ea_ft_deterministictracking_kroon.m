function varargout=ea_ft_deterministictracking_kroon(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Deterministic Fibertracking (Kroon)';
    varargout{2}={'SPM8','SPM12'};
    return
end

%% mask for tracking
directory=[options.root,options.patientname,filesep];
ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0);

%% Coreg options.prefs.prenii_unnormalized to b0 (for label.mat and FTR-Normalization)
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[directory,'c1',options.prefs.prenii_unnormalized,',1'];
                                                   [directory,'c2',options.prefs.prenii_unnormalized,',1'];
                                                   [directory,'c3',options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'rb0';

jobs{1}=matlabbatch;
spm_jobman('run',jobs);
clear matlabbatch jobs;

%% prepare DTI
Vdti=spm_vol([options.root,options.patientname,filesep,options.prefs.dti]);
Xdti=spm_read_vols(Vdti);

bval=load([options.root,options.patientname,filesep,options.prefs.bval]);
bvec=load([options.root,options.patientname,filesep,options.prefs.bvec]);
if size(bval,1)>size(bval,2)
    bval=bval';
end
if size(bvec,1)>size(bvec,2)
    bvec=bvec';
end

for i=1:size(Xdti,4)
   DTIdata(i).VoxelData=single(squeeze(Xdti(:,:,:,i)));
   DTIdata(i).Gradient=bvec(:,i);
   DTIdata(i).Bvalue=bval(:,i);
end

% Constants DTI
parametersDTI=[];
parametersDTI.BackgroundThreshold=150;
parametersDTI.WhiteMatterExtractionThreshold=0.10;
parametersDTI.textdisplay=true;

% Perform DTI calculation
[ADC,FA,VectorF,DifT]=ea_DTI(DTIdata,parametersDTI);

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
parametersFT.Sampling=1;
parametersFT.textdisplay=true;
%parametersFT.FAmin=0.2;
%parametersFT.FAmax=0.8;

% Perform FT
fibers=ea_FT(FA,VectorF,Xmask,parametersFT);
save([options.root,options.patientname,filesep,options.prefs.FTR_unnormalized],'fibers');

%% export .trk copy for trackvis visualization
ea_ftr2trk([directory,options.prefs.FTR_unnormalized],[directory,options.prefs.b0]); % export unnormalized ftr to .trk
disp('Done.');

%% add methods dump:
cits={
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    };
ea_methods(options,['A whole-brain fiber-set was estimated based using a freely-available deterministic diffusion tensor imaging approach as implemented by Dirk Jan-Kroon (https://www.mathworks.com/matlabcentral/fileexchange/21130-dti-and-fiber-tracking).',...
    ' This was done within a white-matter mask that was estimated on the anatomical scan using the Unified Segmentation approach (Ashburner 2005) as implemented in ',spm('ver'),' and linearly co-registered to the b0-weighted series.'],cits);
