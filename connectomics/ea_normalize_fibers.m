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

% check which normalization routine has been used..
% if dartel was used, we need to coregister c2 of b0 and rc2 of anat (since
% deformation fields were estimated for the rc* files and not the native
% anat file.
[options.prefs.b0,options.prefs.prenii_unnormalized,whichnormmethod,template]=ea_checkdartelused(options);

% plot reference volumes
b0=ea_load_nii([directory,options.prefs.b0]);

if vizz
    figure('color','w','name','Fibertrack normalization','numbertitle','off');
    % plot b0
    subplot(1,3,1);
    title('b0 space');
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img>max(b0.img(:))/7));
    plot3(xx(1:10:end),yy(1:10:end),zz(1:10:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    
    % plot anat
    anat=ea_load_nii([directory,options.prefs.prenii_unnormalized]);
    subplot(1,3,2);
    title('anat space');
    [xx,yy,zz]=ind2sub(size(anat.img),find(anat.img>max(anat.img(:))/3));
    plot3(xx(1:1000:end),yy(1:1000:end),zz(1:1000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    
    % plot MNI
    mni=ea_load_nii(template);
    subplot(1,3,3);
    title('MNI space');
    [xx,yy,zz]=ind2sub(size(mni.img),find(mni.img>max(mni.img(:))/3));
    plot3(xx(1:10000:end),yy(1:10000:end),zz(1:10000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
end

% load fibers
[fibers,idx]=ea_loadfibertracts([directory,options.prefs.FTR_unnormalized]);
wfibs=fibers(:,1:3);

% plot unnormalized fibers
if vizz
    try
        thisfib=wfibs(1:100000,:);
    catch
        thisfib=wfibs;
    end
    subplot(1,3,1)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
end

display(sprintf('\nNormalizing fibers...'));

%% Normalize fibers

%% first apply affine transform from b0 to prenii
Vb0=spm_vol([directory,options.prefs.b0]);
Vmprage=spm_vol([directory,options.prefs.prenii_unnormalized]);
x=spm_coreg(Vb0,Vmprage);
affine=Vmprage.mat\spm_matrix(x(:)')*Vb0.mat;
wfibs=[wfibs,ones(size(wfibs,1),1)]*affine';
wfibs=wfibs(:,1:3);

% plot fibers in anat space
if vizz
    try
        thisfib=wfibs(1:100000,1:3);
    catch
        thisfib=wfibs;
    end
    subplot(1,3,2)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
end

%% map from prenii voxelspace to mni millimeter space
switch whichnormmethod
    case {'ea_normalize_spmdartel','ea_normalize_spmnewseg'}
        wfibs = vox2mm_norm(wfibs,[directory,'y_ea_inv_normparams.nii']);
    case ea_getantsnormfuns
        % tranform from anat voxel space to anat mm space
        fibers_mm_mprage=[wfibs,ones(size(wfibs,1),1)]*Vmprage.mat';
        
        % normalize
        wfibs = ea_ants_applytransforms_to_points(directory,fibers_mm_mprage(:,1:3),0);
end

%% map from mni millimeter space to mni voxel space (only needed for trackvis convertion and cleansing fibers).
affine=spm_get_space(template);
wfibsvox=[wfibs,ones(size(wfibs,1),1)]*inv(affine)';
wfibsvox=wfibsvox(:,1:3);

display(sprintf('\nNormalization done.'));

%% cleansing fibers..
if cleanse_fibers % delete anything too far from wm.
    ea_error('Clease fibers not supported at present');
    mnimask=spm_read_vols(spm_vol(template)); % FIX_ME: NEED WM VOLUME OF mni_hires.nii
    mnimask=mnimask>0.01;
    todelete = ~mnimask(sub2ind(size(mnimask),round(wfibsvox(:,1)),round(wfibsvox(:,2)),round(wfibsvox(:,3))));

    wfibs(todelete,:)=[];
    wfibsvox(todelete,:)=[];
end

% plot fibers in MNI space
if vizz
    try
        thisfib=wfibsvox(1:100000,:);
    catch
        thisfib=wfibsvox;
    end
    subplot(1,3,3)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
    drawnow;
end

wfibs=[wfibs,fibers(:,4)];
wfibsvox=[wfibsvox,fibers(:,4)];
[~,ftrbase]=fileparts(options.prefs.FTR_normalized);
ea_savefibertracts([directory,ftrbase,'.mat'],wfibs,idx,'mm');
ea_savefibertracts([directory,ftrbase,'_vox.mat'],wfibsvox,idx,'vox',affine);

%% create normalized trackvis version
try
    display(sprintf('\nExporting normalized fibers to TrackVis...'));
    dnii=ea_load_nii(template);

    specs.origin=[0,0,0];
    specs.dim=size(dnii.img);
    specs.vox=dnii.hdr.dime.pixdim;
    specs.affine=dnii.mat;

    [~,ftrfname]=fileparts(options.prefs.FTR_normalized);
    ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
    disp('Done.');
end


function [useb0,useanat,whichnormmethod,template]=ea_checkdartelused(options)
% check normalization routine used, determine template
useb0=options.prefs.b0;
useanat=options.prefs.prenii_unnormalized;

directory=[options.root,options.patientname,filesep];
whichnormmethod=ea_whichnormmethod(directory);
if isempty(whichnormmethod)
    ea_error('Please run normalization for this subject first.');
end

switch whichnormmethod
    case 'ea_normalize_spmdartel'
        % segment b0.
        if ~exist([directory,'c2',options.prefs.b0],'file');
            disp('Segmenting B0 file for DARTEL import space coregistration...');
            ea_newseg(directory,options.prefs.b0,0,options);
            delete([directory,'c4',options.prefs.b0]);
            delete([directory,'c5',options.prefs.b0]);
            disp('Done.');
        end
        % coreg b0 and anat
        if ~exist([directory,'rc2',options.prefs.prenii_unnormalized],'file');
            ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
            copyfile([directory,options.prefs.prenii_unnormalized],[directory,'k',options.prefs.prenii_unnormalized]);
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,options.prefs.b0]};
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,'k',options.prefs.prenii_unnormalized]};
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[directory,'c1',options.prefs.prenii_unnormalized];
                [directory,'c2',options.prefs.prenii_unnormalized];
                [directory,'c3',options.prefs.prenii_unnormalized]
                };
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
            spm_jobman('run',{matlabbatch});
            clear matlabbatch
        end

        if ~exist([options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6_hires.nii'],'file')
            gunzip([options.earoot,'templates',filesep,'dartel',filesep,'create_mni_darteltemplate',filesep,'dartelmni_6_hires.nii.gz'],...
                   [options.earoot,'templates',filesep,'dartel']);
        end

        useb0=['c2',options.prefs.b0];
        useanat=['rc2',options.prefs.prenii_unnormalized];
        template=[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6_hires.nii,2'];  
	case 'ea_normalize_spmnewseg'
        template=[options.earoot,'templates',filesep,'TPM_Lorio_Draganski.nii,2'];
    otherwise
        template=[options.earoot,'templates',filesep,'mni_hires.nii'];
end


function mm_norm = vox2mm_norm(vox, V)
% returns normalized mm coordinates based on deformation field

if ischar(V)
    V = spm_vol([repmat(V,3,1),[',1,1';',1,2';',1,3']]);
end

vox = double(vox);
mm_norm = [spm_sample_vol(V(1),vox(:,1),vox(:,2),vox(:,3),1);
          spm_sample_vol(V(2),vox(:,1),vox(:,2),vox(:,3),1);
          spm_sample_vol(V(3),vox(:,1),vox(:,2),vox(:,3),1)]';
