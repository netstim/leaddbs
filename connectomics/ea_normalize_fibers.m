function ea_normalize_fibers(options)
% Normalize fiber to MNI space

vizz=1; % turn this value to 1 to visualize fiber normalization (option for debugging only, this will drastically slow down the process).
cleanse_fibers=0; % deletes everything outside the white matter of the template.
directory=[options.root,options.patientname,filesep];

% create unnormalized trackvis version
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
try
	if ~exist([directory,ftrfname,'.trk'],'file')
        display(sprintf('\nExporting unnormalized fibers to TrackVis...'));
        dnii=ea_load_nii([directory,options.prefs.b0]);

        specs.origin=[0,0,0];
        specs.dim=size(dnii.img);
        specs.vox=dnii.hdr.dime.pixdim;
        specs.affine=dnii.mat;

        ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
        disp('Done.');
	end
end

% get transform from b0 to anat and affine matrix of anat
[refb0,refanat,refnorm,b02anat,whichnormmethod]=ea_checktransform(options);

% plot reference volumes
if vizz
    figure('color','w','name',['Fibertrack normalization: ',options.patientname],'numbertitle','off');
    % plot b0
    b0=ea_load_nii(refb0);
    subplot(1,3,1);
    title('b0 space');
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img>max(b0.img(:))/7));
    plot3(xx(1:10:end),yy(1:10:end),zz(1:10:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    
    % plot anat
    anat=ea_load_nii(refanat);
    subplot(1,3,2);
    title('anat space');
    [xx,yy,zz]=ind2sub(size(anat.img),find(anat.img>max(anat.img(:))/3));
    plot3(xx(1:1000:end),yy(1:1000:end),zz(1:1000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    
    % plot MNI
    mni=ea_load_nii(refnorm);
    subplot(1,3,3);
    title('MNI space');
    [xx,yy,zz]=ind2sub(size(mni.img),find(mni.img>max(mni.img(:))/3));
    XYZ_mm=[xx,yy,zz,ones(length(xx),1)]*mni.mat';
    plot3(XYZ_mm(1:10000:end,1),XYZ_mm(1:10000:end,2),XYZ_mm(1:10000:end,3),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
end

% load fibers
[fibers,idx]=ea_loadfibertracts([directory,options.prefs.FTR_unnormalized]);

% plot unnormalized fibers
if vizz
    try
        thisfib=fibers(1:100000,:);
    catch
        thisfib=fibers;
    end
    subplot(1,3,1)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
end

display(sprintf('\nNormalizing fibers...'));

%% Normalize fibers

%% map from b0 voxel space to anat voxel space
wfibsvox_anat=[fibers(:,1:3),ones(size(fibers,1),1)]*b02anat';
wfibsvox_anat=wfibsvox_anat(:,1:3);

%% map from anat voxel space to anat mm space
anataffine=spm_get_space(refanat);
wfibsmm_anat=[wfibsvox_anat,ones(size(wfibsvox_anat,1),1)]*anataffine';
wfibsmm_anat=wfibsmm_anat(:,1:3);

% plot fibers in anat space
if vizz
    try
        thisfib=wfibsvox_anat(1:100000,:);
    catch
        thisfib=wfibsvox_anat;
    end
    subplot(1,3,2)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
end

%% map from anat voxel space to mni mm space
display(sprintf('\nPoints normalization...'));
wfibsmm_mni = ea_map_coords(wfibsvox_anat',refanat,[directory,'y_ea_inv_normparams.nii'])';

%% map from mni mm space to mni voxel space
mniaffine=spm_get_space(refnorm);
wfibsvox_mni=[wfibsmm_mni,ones(size(wfibsmm_mni,1),1)]*inv(mniaffine)';
wfibsvox_mni=wfibsvox_mni(:,1:3);

display(sprintf('\nNormalization done.'));

%% cleansing fibers..
if cleanse_fibers % delete anything too far from wm.
    ea_error('Clease fibers not supported at present');
    mnimask=spm_read_vols(spm_vol(refnorm)); % FIX_ME: NEED WM VOLUME OF mni_hires.nii
    mnimask=mnimask>0.01;
    todelete = ~mnimask(sub2ind(size(mnimask),round(wfibsvox_mni(:,1)),round(wfibsvox_mni(:,2)),round(wfibsvox_mni(:,3))));

    wfibsmm_mni(todelete,:)=[];
    wfibsvox_mni(todelete,:)=[];
end

% plot fibers in MNI space
if vizz
    try
        thisfib=wfibsmm_mni(1:100000,:);
    catch
        thisfib=wfibsmm_mni;
    end
    subplot(1,3,3)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
    drawnow;
end

%% export fibers
wfibsmm_mni=[wfibsmm_mni,fibers(:,4)];
wfibsvox_mni=[wfibsvox_mni,fibers(:,4)];
[~,ftrbase]=fileparts(options.prefs.FTR_normalized);
if ~exist([directory,'connectomes',filesep,'dMRI'],'file')
    mkdir([directory,'connectomes',filesep,'dMRI']);
end
ea_savefibertracts([directory,'connectomes',filesep,'dMRI',filesep,ftrbase,'.mat'],wfibsmm_mni,idx,'mm');
ea_savefibertracts([directory,'connectomes',filesep,'dMRI',filesep,ftrbase,'_vox.mat'],wfibsvox_mni,idx,'vox',mniaffine);

%% create normalized trackvis version
try
    display(sprintf('\nExporting normalized fibers to TrackVis...'));
    dnii=ea_load_nii(refnorm);

    specs.origin=[0,0,0];
    specs.dim=size(dnii.img);
    specs.vox=dnii.hdr.dime.pixdim;
    specs.affine=dnii.mat;

    [~,ftrfname]=fileparts(options.prefs.FTR_normalized);
    ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
    disp('Done.');
end


function [refb0,refanat,refnorm,b02anat,whichnormmethod]=ea_checktransform(options)
directory=[options.root,options.patientname,filesep];
% % segment b0.
% if ~exist([directory,'c2',options.prefs.b0],'file');
%     disp('Segmenting b0...');
%     ea_newseg(directory,options.prefs.b0,0,options);
%     disp('Done.');
% end

% % segment anat.
% if ~exist([directory,'c2',options.prefs.prenii_unnormalized],'file');
%     disp('Segmenting anat...');
%     ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
%     disp('Done.');
% end

% check normalization routine used, determine template
[whichnormmethod,refnorm]=ea_whichnormmethod(directory);
if isempty(whichnormmethod)
    ea_error('Please run normalization for this subject first.');
end

% determine the refimage for b0 and anat space visualization
refb0=[directory,options.prefs.b0];
refanat=[directory,options.prefs.prenii_unnormalized];

% determin the template for fiber normalization and visualization
if ismember(whichnormmethod,{'ea_normalize_spmdartel','ea_normalize_spmnewseg'})
	refnorm=[refnorm,',2'];  
end

% unzip spm dartel template
if strcmp(whichnormmethod, 'ea_normalize_spmdartel')
    if ~exist([options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6_hires.nii'],'file')
        gunzip([options.earoot,'templates',filesep,'dartel',filesep,'create_mni_darteltemplate',filesep,'dartelmni_6_hires.nii.gz'],...
               [options.earoot,'templates',filesep,'dartel']);
    end
end

% generate b0 to anat tranformation
Vb0=spm_vol([directory,options.prefs.b0]);
Vanat=spm_vol([directory,options.prefs.prenii_unnormalized]);

switch options.coregmr.method
    case 2 % ANTs
        ea_ants([directory,options.prefs.prenii_unnormalized],[directory,options.prefs.b0],[directory,'r',options.prefs.b0],1);
        load([directory,'ct2anat1.mat']);
        delete([directory,'ct2anat1.mat']);
        b02anat=eye(4);
        b02anat(1:3,1)=AffineTransform_float_3_3(1:3);
        b02anat(1:3,2)=AffineTransform_float_3_3(4:6);
        b02anat(1:3,3)=AffineTransform_float_3_3(7:9);
        b02anat(1:3,4)=AffineTransform_float_3_3(10:12);

        % the following is only empirically determined for now. This could
        % be wrong in some cases.

        b02anat([7,10,15])=b02anat([7,10,15])*-1;
        b02anat=Vanat.mat\b02anat*Vb0.mat;
    otherwise % default use SPM, BRAINSFit not supported here
        display('Register b0 to anat...');
        x=spm_coreg(Vb0,Vanat);
        b02anat=Vanat.mat\spm_matrix(x(:)')*Vb0.mat;
end
save([directory,'b02anat.mat'],'b02anat');
