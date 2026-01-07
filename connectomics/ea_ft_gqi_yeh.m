function varargout=ea_ft_gqi_yeh(options)
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

if ischar(options) % return name of method.
    varargout{1}='Generalized Q-Sampling (Yeh et al. 2010)';
    varargout{2}={'SPM8','SPM12'};
    return
end

directory=[options.root,options.patientname,filesep];

[redo, options]=ea_prepare_dti(options);

vizz=0;

% build white matter mask
% Check both classic location (root) and BIDS location (preprocessing/dwi/)
maskClassic = [directory,'ttrackingmask.nii'];
maskBIDS = fullfile(directory, 'preprocessing', 'dwi', 'ttrackingmask.nii');
if (~isfile(maskClassic) && ~isfile(maskBIDS)) || redo || ...
        (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    ea_gentrackingmask_brainmask(options,1)
end

basedir = [options.earoot, 'ext_libs',filesep,'dsi_studio',filesep];
DSISTUDIO = ea_getExec([basedir, 'dsi_studio'], 'escapePath', 1);

% build .fz file
[~,ftrbase]=fileparts(options.prefs.FTR_unnormalized);
% BIDS: Check for .fz file in connectomics/dMRI/
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
fibFile = fullfile(connectomicsDir, [ftrbase, '.fz']);
if ~isfile(fibFile) || redo || ...
        (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    disp('Estimating ODF / preparing GQI...');

    ea_prepare_fib_gqi(DSISTUDIO,options);

    disp('Done.');
else
    disp('.fz file found, no need to rebuild.');
end

% BIDS: Use paths from connectomics/dMRI/ for FIB and output
trkcmd=[DSISTUDIO,' --action=trk',...
    ' --method=0',...
    ' --source=',ea_path_helper(fullfile(connectomicsDir, [ftrbase,'.fz'])),...
    ' --seed=',ea_path_helper([directory,'ttrackingmask.nii']),...
    ' --fiber_count=', num2str(options.lc.struc.ft.dsistudio.fiber_count),...
    ' --output=',ea_path_helper(fullfile(connectomicsDir, [ftrbase,'.mat'])),...
    ' --dt_threshold=0.2',...
    ' --max_length=300.0',...
    ' --min_length=10.0',...
    ' --random_seed=0',...
    ' --smoothing=0.2',...
    ' --step_size=0.5',...
    ' --turning_angle=75'];

err=ea_runcmd(trkcmd);
if err
    ea_error(['Fibertracking with dsi_studio failed (error code=',num2str(err),').']);
end

ea_delete([directory,'ttrackingmask.txt']);
% now store tract in lead-dbs format
disp('Converting fibers...');
fibinfo = load(fullfile(connectomicsDir, [ftrbase,'.mat']));

% DEBUG: Check DSI Studio raw output
fprintf('\n=== DEBUG: DSI Studio RAW OUTPUT ===\n');
fprintf('Fields: %s\n', strjoin(fieldnames(fibinfo), ', '));
fprintf('tracts size: [%d, %d]\n', size(fibinfo.tracts, 1), size(fibinfo.tracts, 2));
fprintf('tracts range (DSI format, 0-based):\n');
fprintf('  X: [%.2f, %.2f]\n', min(fibinfo.tracts(:,1)), max(fibinfo.tracts(:,1)));
fprintf('  Y: [%.2f, %.2f]\n', min(fibinfo.tracts(:,2)), max(fibinfo.tracts(:,2)));
fprintf('  Z: [%.2f, %.2f]\n', min(fibinfo.tracts(:,3)), max(fibinfo.tracts(:,3)));

fibers = fibinfo.tracts';
idx = double(fibinfo.length)';
clear fibinfo
b0 = spm_vol([directory,options.prefs.b0]);

fprintf('\nTransposed fibers size: [%d, %d]\n', size(fibers, 1), size(fibers, 2));
fprintf('b0.dim = [%d, %d, %d]\n', b0.dim(1), b0.dim(2), b0.dim(3));

fprintf('DSI output range (0-based): X=[%.2f, %.2f], Y=[%.2f, %.2f], Z=[%.2f, %.2f]\n', ...
    min(fibers(:,1)), max(fibers(:,1)), min(fibers(:,2)), max(fibers(:,2)), min(fibers(:,3)), max(fibers(:,3)));

% Default orientation in DSI-Studio and TrackVis is LPS. Flip the
% coordinates to make the orientation in the MAT file inline with b0 image.
if b0.mat(1)>0  % 'R' is positive x-axis
    % flip x
    disp('Flip positive X-axis to R...');
    fibers(:,1) = b0.dim(1)-1-fibers(:,1);
end
if b0.mat(6)>0  % 'A' is positive y-axis
    %flip y
    disp('Flip positive Y-axis to A...');
    fibers(:,2) = b0.dim(2)-1-fibers(:,2);
end
if b0.mat(11)<0  % 'I' is positive z-axis
    %flip z
    disp('Flip positive Z-axis to I...');
    fibers(:,3) = b0.dim(3)-1-fibers(:,3);
end

% Change ZERO-BASED indexing to ONE-BASED indexing.
fibers(:,1:3) = fibers(:,1:3) + 1;

fprintf('After flip and 1-based conversion:\n');
fprintf('  X: [%.2f, %.2f] (expected: [1, %d])\n', min(fibers(:,1)), max(fibers(:,1)), b0.dim(1));
fprintf('  Y: [%.2f, %.2f] (expected: [1, %d])\n', min(fibers(:,2)), max(fibers(:,2)), b0.dim(2));
fprintf('  Z: [%.2f, %.2f] (expected: [1, %d])\n', min(fibers(:,3)), max(fibers(:,3)), b0.dim(3));

if vizz
    figure
    thresh=700; % set to a good grey value.
    plot3(fibers(:,1),fibers(:,2),fibers(:,3),'r.')
    hold on
    b0=ea_load_nii([directory,options.prefs.b0]);
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img(:)>thresh));
    plot3(xx,yy,zz,'g.')
end

% BIDS FIX: Save fibers as [N x 3] matrix, not [N x 4]
% The fiber index is stored separately in 'idx'
ftr.ea_fibformat = '1.0';
ftr.fourindex = 0;  % No 4th column
ftr.fibers = fibers(:, 1:3);  % Only X, Y, Z
ftr.idx = idx;
ftr.voxmm = 'vox';
ftr.mat = b0.mat;

% BIDS: Save fibers in connectomics/dMRI/
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
if ~exist(connectomicsDir, 'dir')
    mkdir(connectomicsDir);
end

disp('Saving fibers...');
save(fullfile(connectomicsDir, [ftrbase, '.mat']), '-struct', 'ftr', '-v7.3');
disp('Done.');

fprintf('\nGenerating trk in b0 space...\n');
ea_ftr2trk(fullfile(connectomicsDir, [ftrbase, '.mat']), [directory, options.prefs.b0])


function ea_prepare_fib_gqi(DSISTUDIO,options)
directory=[options.root,options.patientname,filesep];
[~,ftrbase]=fileparts(options.prefs.FTR_unnormalized);

% Create SRC file
srcFile = [directory,'dti.sz'];
ea_delete(srcFile);
cmd = [DSISTUDIO,' --action=src --source=',ea_path_helper([directory,options.prefs.dti]),...
       ' --bval=',ea_path_helper([directory,options.prefs.bval]),...
       ' --bvec=',ea_path_helper([directory,options.prefs.bvec]),...
       ' --output=',ea_path_helper(srcFile),...
       ' --sort_b_table=0'];

err = ea_runcmd(cmd);

if err || ~isfile(srcFile)
    ea_error('DSI studio failed to generate the SRC file!', simpleStack=1);
end

% Create FIB file (BIDS: save in connectomics/dMRI/)
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
if ~exist(connectomicsDir, 'dir')
    mkdir(connectomicsDir);
end
fibFile = fullfile(connectomicsDir, [ftrbase, '.fz']);
cmd = [DSISTUDIO,' --action=rec --source=',ea_path_helper(srcFile),...
       ' --dti_no_high_b=1',...
       ' --mask=',ea_path_helper([directory,'ttrackingmask.nii']),...
       ' --method=4',...
       ' --param0=1.25',...
       ' --check_btable=0',...
       ' --make_isotropic=0',...
       ' --record_odf=0',...
       ' --other_output=fa,rd,iso,rdi',...
       ' --output=',ea_path_helper(fibFile)];

err = ea_runcmd(cmd);
ea_delete(srcFile);

if err || ~isfile(fibFile)
    ea_error('DSI studio failed to generate the FIB file!', simpleStack=1);
end

%% add methods dump:
cits={
    'Yeh, F.-C., Wedeen, V. J., & Tseng, W.-Y. I. (2010). Generalized q-sampling imaging. IEEE Transactions on Medical Imaging, 29(9), 1626?1635. http://doi.org/10.1109/TMI.2010.2045126'
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    };
ea_methods(options,['A whole-brain fiber-set was estimated based using the Generalized q-sampling imaging (GQI) approach (Yeh 2010) as implemented in DSI-Studio (http://dsi-studio.labsolver.org).',...
    ' GQI is a model-free method that calculates the orientational distribution of the density of diffusing water.',...
    ' Fibers were sampled within a white-matter mask that was estimated using the anatomical acquisition by applying the Unified Segmentation approach (Ashburner 2005) as implemented in ',spm('ver'),'. This mask was linearly co-registered to the b0-weighted series.'],cits);
