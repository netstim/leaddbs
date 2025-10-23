function ea_normalize_fibers(options)
% Normalize fiber to MNI space

vizz=1; % turn this value to 1 to visualize fiber normalization (option for debugging only, this will drastically slow down the process).
cleanse_fibers=0; % deletes everything outside the white matter of the template.
directory=[options.root,options.patientname,filesep];

% create unnormalized trackvis version
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
connectomicsDir = fullfile(directory, 'connectomics', 'dMRI');
try
	if ~exist(fullfile(connectomicsDir, [ftrfname, '.trk']), 'file')
        fprintf('\nExporting unnormalized fibers to TrackVis...\n');
        ea_ftr2trk(fullfile(connectomicsDir, [ftrfname, '.mat']), [directory, options.prefs.b0]);
        disp('Done.');
	end
end

% BIDS FIX: Check if normalization exists - if not, run SPM12 normalization
normFile = fullfile(directory, 'y_ea_inv_normparams.nii');
if ~exist(normFile, 'file')
    fprintf('\nNo normalization found for fiber normalization. Running SPM12 normalization...\n');
    
    % Get anatomical image
    anatFile = fullfile(directory, options.prefs.prenii_unnormalized);
    
    % Run SPM12 normalization (quick estimate)
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {[anatFile, ',1']};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {[anatFile, ',1']};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {[spm('Dir'), '/tpm/TPM.nii']};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    
    spm_jobman('run', {matlabbatch});
    clear matlabbatch
    
    % Rename output to Lead-DBS convention
    [anatDir, anatName] = fileparts(anatFile);
    spmNormFile = fullfile(anatDir, ['y_', anatName, '.nii']);
    if exist(spmNormFile, 'file')
        movefile(spmNormFile, normFile);
    end
    
    fprintf('Normalization complete.\n\n');
end

% get transform from b0 to anat and affine matrix of anat
[refb0,refanat,refnorm,whichnormmethod]=ea_checktransform(options);

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
[~, mov] = fileparts(options.prefs.b0);
[~, fix] = fileparts(options.prefs.prenii_unnormalized);

% For BIDS: search in coregistration/transformations and preprocessing dirs
if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    % BIDS FIX: Use dir() instead of ea_regexpdir for ANTs files
    searchDirs = {
        fullfile(directory, 'coregistration', 'transformations')
        fullfile(directory, 'preprocessing', 'anat')
        directory
    };
    
    transform = {};
    for iDir = 1:length(searchDirs)
        if exist(searchDirs{iDir}, 'dir')
            % Search for ANTs composite files
            pattern = [mov, '2', fix, '*Composite.nii.gz'];
            files = dir(fullfile(searchDirs{iDir}, pattern));
            if ~isempty(files)
                % Return both forward and inverse if they exist
                for iFile = 1:length(files)
                    transform{end+1} = fullfile(searchDirs{iDir}, files(iFile).name);
                end
                break;
            end
        end
    end
else
    % Extract just the method name without version/citation info
    coregmethod = strrep(options.coregmr.method, 'Hybrid SPM & ', '');
    % Remove everything after and including the first space or parenthesis
    coregmethod = regexp(coregmethod, '^[^\s\(]+', 'match', 'once');
    if isempty(coregmethod)
        coregmethod = 'spm'; % fallback
    end
    options.coregmr.method = coregmethod;
    
    % BIDS FIX: Use simpler search with dir() instead of ea_regexpdir
    % Search for b0→anat transformation file
    [~, b0Name] = ea_niifileparts(options.prefs.b0);
    [~, anatName] = ea_niifileparts(options.prefs.prenii_unnormalized);
    
    fprintf('DEBUG ea_normalize_fibers: Searching for transformation...\n');
    fprintf('  b0Name: %s\n', b0Name);
    fprintf('  anatName: %s\n', anatName);
    fprintf('  coregmethod: %s\n', lower(coregmethod));
    
    % Try preprocessing/anat first
    searchDirs = {fullfile(directory, 'preprocessing', 'anat'), directory};
    transform = {};
    
    for iDir = 1:length(searchDirs)
        if exist(searchDirs{iDir}, 'dir')
            % Search for files matching pattern: b0Name2anatName_method.mat
            pattern = [b0Name, '2', anatName, '_', lower(coregmethod), '*.mat'];
            fprintf('  Searching in: %s\n', searchDirs{iDir});
            fprintf('  Pattern: %s\n', pattern);
            files = dir(fullfile(searchDirs{iDir}, pattern));
            fprintf('  Found %d files\n', length(files));
            if ~isempty(files)
                transform = {fullfile(searchDirs{iDir}, files(end).name)};
                fprintf('  Using: %s\n', transform{1});
                break;
            end
        end
    end
end

if isempty(transform)
    fprintf('ERROR: No transformation found!\n');
    fprintf('Available files in preprocessing/anat:\n');
    if exist(fullfile(directory, 'preprocessing', 'anat'), 'dir')
        allFiles = dir(fullfile(directory, 'preprocessing', 'anat', '*.mat'));
        for i = 1:length(allFiles)
            fprintf('  %s\n', allFiles(i).name);
        end
    end
    ea_error('Coregistration transformation not found. Please ensure coregistration between b0 and anatomical has been completed before running fiber normalization.');
end

if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    % For ANTs: use the transform we found
    if ~isempty(transform)
        xfmPath = transform{1};
        % Find the InverseComposite version
        xfmPath = regexprep(xfmPath, 'Composite', 'InverseComposite');
        if ~exist(xfmPath, 'file')
            % If InverseComposite doesn't exist, construct the path
            xfmDir = fullfile(directory, 'coregistration', 'transformations');
            xfmPath = fullfile(xfmDir, [mov, '2', fix, 'InverseComposite.nii.gz']);
        end
    else
        % Fallback: construct expected path
        xfmDir = fullfile(directory, 'coregistration', 'transformations');
        xfmPath = fullfile(xfmDir, [mov, '2', fix, 'InverseComposite.nii.gz']);
    end
    [~, wfibsvox_anat] = ea_map_coords(fibers(:,1:3)', ...
                                       refb0, ...
                                       xfmPath, ...
                                       refanat, 'ANTs');
else
    % For SPM/FSL: use .mat file
    if ~isempty(transform)
        xfmPath = transform{1};
    else
        % Fallback: construct expected path in preprocessing/anat
        xfmPath = fullfile(directory, 'preprocessing', 'anat', [mov, '2', fix, '.mat']);
    end
    [~, wfibsvox_anat] = ea_map_coords(fibers(:,1:3)', ...
                                       refb0, ...
                                       xfmPath, ...
                                       refanat, ...
                                       options.coregmr.method);
end

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

% Get BIDS-compliant transformation files
transformfiles = ea_gettransformfiles(options);

% Use inverse transform (from native to MNI)
% Note: ea_map_coords will add the appropriate extension if needed
if exist(transformfiles.inverse, 'file')
    inverseTransform = transformfiles.inverse;
else
    % Fallback to classic naming (ea_map_coords will add extension)
    inverseTransform = [directory,'inverseTransform'];
end

% Determine transformation method from whichnormmethod
if contains(whichnormmethod, 'ANTs', 'IgnoreCase', true)
    transformmethod = 'ANTs';
elseif contains(whichnormmethod, 'SPM', 'IgnoreCase', true)
    transformmethod = 'SPM';
elseif contains(whichnormmethod, 'FNIRT', 'IgnoreCase', true)
    transformmethod = 'FSL';
else
    transformmethod = 'ANTs'; % default
end

[wfibsmm_mni, wfibsvox_mni] = ea_map_coords(wfibsvox_anat', ...
                                            refanat, ...
                                            inverseTransform, ...
                                            refnorm, ...
                                            transformmethod);

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

% BIDS FIX: If no method detected but y_ea_inv_normparams.nii exists, assume SPM12
if isempty(whichnormmethod)
    normFile = fullfile(directory, 'y_ea_inv_normparams.nii');
    if exist(normFile, 'file')
        fprintf('Normalization file found, assuming SPM12 method...\n');
        whichnormmethod = 'SPM12';
        % Use default template
        spacedef = ea_getspacedef;
        refnorm = [ea_space, spacedef.templates{1}, '.nii'];
    else
        ea_error('Please run normalization for this subject first.');
    end
end

% determine the refimage for b0 and anat space visualization
refb0=[directory,options.prefs.b0];
refanat=[directory,options.prefs.prenii_unnormalized];

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
