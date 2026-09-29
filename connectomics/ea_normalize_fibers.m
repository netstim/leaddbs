function ea_normalize_fibers(options)
% Normalize fiber to MNI space

vizz=1; % turn this value to 1 to visualize fiber normalization (option for debugging only, this will drastically slow down the process).
cleanse_fibers=0; % deletes everything outside the white matter of the template.
directory=[options.root,options.patientname,filesep];

record = ea_lc_coreg_record(options);
nativeB0 = record.moving;
options.prefs.b0 = erase(record.moving, directory);
options.prefs.prenii_unnormalized = erase(record.anchor, directory);

% create unnormalized trackvis version
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
try
	if ~exist(fullfile(connectomicsDir, [ftrfname, '.trk']), 'file')
        fprintf('\nExporting unnormalized fibers to TrackVis...\n');
        ea_ftr2trk(fullfile(connectomicsDir, [ftrfname, '.mat']), nativeB0);
        disp('Done.');
	end
end

% Require the active normalization protocol. Never estimate a second warp here.
transformfiles = ea_gettransformfiles(options);
if ~isfile(transformfiles.forward) || ~isfile(transformfiles.inverse)
    error('LeadDBS:MissingNormalization', 'Run anatomical normalization before normalizing fibers.');
end

% get transform from b0 to anat and affine matrix of anat
[refb0,refanat,refnorm,whichnormmethod]=ea_checktransform(options);
refb0 = record.moving;
refanat = record.anchor;

% plot reference volumes
if vizz
    figure('color','w','name',['Fibertrack normalization: ',options.patientname],'numbertitle','off');
    % plot b0, voxel space
    b0=ea_load_nii(refb0);
    subplot(1,3,1);
    title('b0 space');
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img>max(b0.img(:))/7));
    plot3(xx(1:10:end),yy(1:10:end),zz(1:10:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    title('b0 space');
    hold on

    % plot anat, voxel space
    anat=ea_load_nii(refanat);
    subplot(1,3,2);
    title('anat space');
    [xx,yy,zz]=ind2sub(size(anat.img),find(anat.img>max(anat.img(:))/3));
    plot3(xx(1:1000:end),yy(1:1000:end),zz(1:1000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    title('anat space');
    hold on

    % plot MNI, world space
    mni=ea_load_nii(refnorm);
    subplot(1,3,3);
    title('MNI space');
    [xx,yy,zz]=ind2sub(size(mni.img),find(mni.img>max(mni.img(:))/3));
    XYZ_mm=[xx,yy,zz,ones(length(xx),1)]*mni.mat';
    plot3(XYZ_mm(1:10000:end,1),XYZ_mm(1:10000:end,2),XYZ_mm(1:10000:end,3),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    title('MNI space');
    hold on
end

% load fibers
% BIDS: Load from connectomics/dMRI/
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
ftrPath = fullfile(connectomicsDir, [ftrfname, '.mat']);
if ~exist(ftrPath, 'file')
    % Fallback: try classic root location
    ftrPath = [directory, options.prefs.FTR_unnormalized];
end
[fibers,idx]=ea_loadfibertracts(ftrPath);

% plot unnormalized fibers
maxvisfiber = 100000;
if vizz
    if size(fibers,1) > maxvisfiber
        thisfib=fibers(1:maxvisfiber,:);
    else
        thisfib=fibers;
    end
    subplot(1,3,1)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
end

fprintf('\nNormalizing fibers...\n');

%% Normalize fibers

%% map from b0 voxel space to anat mm and voxel space
fprintf('\nMapping from b0 to anat...\n');

% Use the recorded transform and its source grid.
fprintf('Using recorded B0-to-anatomy transform: %s (%s)\n', record.forward, record.method);
if record.nonlinear
    mappedMM = ea_map_coords(fibers(:,1:3)', record.moving, ...
        record.inverse, record.fixed, 'ANTs', 0);
else
    mappedMM = ea_map_coords(fibers(:,1:3)', record.moving, ...
        record.forward, record.fixed, record.method);
end
% Alternate fixed anatomy is already in anchorNative world space.
wfibsvox_anat = ea_get_affine(refanat) \ [mappedMM; ones(1, size(mappedMM,2))];
wfibsvox_anat = wfibsvox_anat(1:3,:);

wfibsvox_anat = wfibsvox_anat';

% DEBUG: Check transformation results
fprintf('DEBUG: Original fibers size: %d x %d\n', size(fibers));
fprintf('DEBUG: Transformed anat fibers size: %d x %d\n', size(wfibsvox_anat));
fprintf('DEBUG: Anat fibers range: X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f]\n', ...
    min(wfibsvox_anat(:,1)), max(wfibsvox_anat(:,1)), ...
    min(wfibsvox_anat(:,2)), max(wfibsvox_anat(:,2)), ...
    min(wfibsvox_anat(:,3)), max(wfibsvox_anat(:,3)));

% BIDS: Save intermediate anat-space fibers in connectomics/dMRI/
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
if ~exist(connectomicsDir, 'dir')
    mkdir(connectomicsDir);
end

% BIDS FIX: Fibers are now [N x 3], no 4th column
ea_savefibertracts(fullfile(connectomicsDir, [ftrfname, '_anat.mat']), wfibsvox_anat, idx, 'vox', refanat);
fprintf('\nGenerating trk in anat space...\n');
ea_ftr2trk(fullfile(connectomicsDir, [ftrfname, '_anat.mat']), refanat);

% plot fibers in anat space
if vizz
    if size(wfibsvox_anat,1) > maxvisfiber
        thisfib=wfibsvox_anat(1:maxvisfiber,:);
    else
        thisfib=wfibsvox_anat;
    end
    subplot(1,3,2)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
end

%% map from anat voxel space to mni mm and voxel space
fprintf('\nMapping from anat to mni...\n');

[wfibsmm_mni, wfibsvox_mni] = ea_lc_anat_to_mni(options, ...
    wfibsvox_anat', refanat, refnorm);

wfibsmm_mni = wfibsmm_mni';
wfibsvox_mni = wfibsvox_mni';

% DEBUG: Check MNI transformation results
fprintf('DEBUG: MNI mm fibers size: %d x %d\n', size(wfibsmm_mni));
fprintf('DEBUG: MNI mm fibers range: X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f]\n', ...
    min(wfibsmm_mni(:,1)), max(wfibsmm_mni(:,1)), ...
    min(wfibsmm_mni(:,2)), max(wfibsmm_mni(:,2)), ...
    min(wfibsmm_mni(:,3)), max(wfibsmm_mni(:,3)));
fprintf('DEBUG: MNI vox fibers size: %d x %d\n', size(wfibsvox_mni));

fprintf('\nNormalization done.\n');

%% cleansing fibers..
if cleanse_fibers % delete anything too far from wm.
    ea_error('Clease fibers not supported at present');
    mnimask=spm_read_vols(spm_vol(refnorm)); % FIX_ME: NEED WM VOLUME OF mni_hires_t2.nii
    mnimask=mnimask>0.01;
    todelete = ~mnimask(sub2ind(size(mnimask),round(wfibsvox_mni(:,1)),round(wfibsvox_mni(:,2)),round(wfibsvox_mni(:,3))));

    wfibsmm_mni(todelete,:)=[];
    wfibsvox_mni(todelete,:)=[];
end

% plot fibers in MNI space
if vizz
    if size(wfibsmm_mni,1) > maxvisfiber
        thisfib=wfibsmm_mni(1:maxvisfiber,:);
    else
        thisfib=wfibsmm_mni;
    end
    subplot(1,3,3)
    plot3(thisfib(:,1),thisfib(:,2),thisfib(:,3),'.','color',[0.1707    0.2919    0.7792]);
    drawnow;
end

%% export fibers
[~,ftrbase]=fileparts(options.prefs.FTR_normalized);
if ~exist([directory,'connectomics',filesep,'dMRI'],'file')
    mkdir([directory,'connectomics',filesep,'dMRI']);
end
% BIDS FIX: Fibers are now [N x 3], no 4th column
ea_savefibertracts([directory,'connectomics',filesep,'dMRI',filesep,ftrbase,'.mat'], wfibsmm_mni, idx, 'mm');
ea_savefibertracts([directory,'connectomics',filesep,'dMRI',filesep,ftrbase,'_vox.mat'], wfibsvox_mni, idx, 'vox', refnorm);

%% create normalized trackvis version
fprintf('\nExporting normalized fibers to TrackVis...\n');

[~,ftrfname]=fileparts(options.prefs.FTR_normalized);
ea_ftr2trk([directory,'connectomics',filesep,'dMRI',filesep,ftrfname]); % export normalized ftr to .trk
disp('Done.');

%% add methods dump:
cits={
    'Horn, A., Ostwald, D., Reisert, M., & Blankenburg, F. (2014). The structural-functional connectome and the default mode network of the human brain. NeuroImage, 102 Pt 1, 142-151. http://doi.org/10.1016/j.neuroimage.2013.09.069'
    'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127-135. http://doi.org/10.1016/j.neuroimage.2014.12.002'
    'Horn, A., & Blankenburg, F. (2016). Toward a standardized structural-functional group connectome in MNI space. NeuroImage, 124(Pt A), 310-322. http://doi.org/10.1016/j.neuroimage.2015.08.048'
    };
ea_methods(options,['The whole-brain fiber set was normalized into standard-stereotactic space following the approach described in (Horn 2014, Horn 2016) as ',...
    ' implemented in Lead-DBS software (Horn 2015; www.lead-dbs.org).'],...
    cits);


function [refb0,refanat,refnorm,whichnormmethod]=ea_checktransform(options)
directory=[options.root,options.patientname,filesep];

% check normalization routine used, determine template
[whichnormmethod,refnorm]=ea_whichnormmethod(directory);

if isempty(whichnormmethod)
    error('LeadDBS:MissingNormalization', 'The active normalization method must be recorded.');
end

% determine the refimage for b0 and anat space visualization
% Primary definition (classic Lead-DBS behaviour)
refb0 = [directory, options.prefs.b0];
refanat = [directory, options.prefs.prenii_unnormalized];

% determine the template for fiber normalization and visualization
if ismember(whichnormmethod,{'ea_normalize_spmshoot','ea_normalize_spmdartel','ea_normalize_spmnewseg'})
	refnorm=[refnorm,',2'];
end

% BIDS fix: If template doesn't exist, use brainmask as fallback
if ~exist(refnorm, 'file')
    % Try .gz version
    if exist([refnorm, '.gz'], 'file')
        refnorm = [refnorm, '.gz'];
    else
        % Fallback to brainmask (same MNI space, just different modality)
        % Extract space name from refnorm path
        [refnorm_dir, refnorm_name] = fileparts(refnorm);
        brainmask = fullfile(refnorm_dir, 'brainmask.nii.gz');
        if exist(brainmask, 'file')
            fprintf('Template %s not found, using brainmask as reference for coordinate transformation.\n', refnorm_name);
            refnorm = brainmask;
        else
            ea_error('Template file not found: %s. Please download Lead-DBS templates.', refnorm);
        end
    end
end
