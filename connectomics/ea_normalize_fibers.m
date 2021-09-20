function ea_normalize_fibers(options)
% Normalize fiber to MNI space

vizz=1; % turn this value to 1 to visualize fiber normalization (option for debugging only, this will drastically slow down the process).
cleanse_fibers=0; % deletes everything outside the white matter of the template.
directory=[options.root,options.patientname,filesep];

% create unnormalized trackvis version
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
try
	if ~exist([directory,ftrfname,'.trk'],'file')
        fprintf('\nExporting unnormalized fibers to TrackVis...\n');
        ea_ftr2trk([directory,ftrfname,'.mat'],[directory,options.prefs.b0]);
        disp('Done.');
	end
end

% get transform from b0 to anat and affine matrix of anat
[refb0,refanat,refnorm]=ea_checktransform(options);

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
[fibers,idx]=ea_loadfibertracts([directory,options.prefs.FTR_unnormalized]);

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
if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    xfm = [mov, '2', fix, '(Inverse)?Composite\.nii\.gz$'];
else
    coregmethod = strrep(options.coregmr.method, 'Hybrid SPM & ', '');
    options.coregmr.method = coregmethod;
    xfm = [mov, '2', fix, '_', lower(coregmethod), '\d*\.(mat|h5)$'];
end
transform = ea_regexpdir(directory, xfm, 0);

if numel(transform) == 0
    warning('Specified transformation not found! Running coregistration now!');
    if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
        ea_ants_nonlinear_coreg([directory,options.prefs.prenii_unnormalized],...
            [directory,options.prefs.b0],...
            [directory,ea_stripext(options.prefs.b0), '2', options.prefs.prenii_unnormalized]);
        ea_delete([directory,ea_stripext(options.prefs.b0), '2', options.prefs.prenii_unnormalized]);
    else
        ea_backuprestore(refb0);
        ea_coregimages(options,refb0,refanat,[options.root,options.patientname,filesep,'tmp.nii'],{},1);
        ea_delete([options.root,options.patientname,filesep,'tmp.nii']);
    end
end

if strcmp(options.coregmr.method, 'ANTs') && options.coregb0.addSyN
    [~, wfibsvox_anat] = ea_map_coords(fibers(:,1:3)', ...
                                       refb0, ...
                                       [directory,ea_stripext(options.prefs.b0), '2', ea_stripext(options.prefs.prenii_unnormalized), 'InverseComposite.nii.gz'], ...
                                       refanat, 'ANTs');
else
    [~, wfibsvox_anat] = ea_map_coords(fibers(:,1:3)', ...
                                       refb0, ...
                                       [directory, mov, '2', fix, '.mat'], ...
                                       refanat, ...
                                       options.coregmr.method);
end

wfibsvox_anat = wfibsvox_anat';
ea_savefibertracts([directory,ftrfname,'_anat.mat'],[wfibsvox_anat,fibers(:,4)],idx,'vox',refanat);
fprintf('\nGenerating trk in anat space...\n');
ea_ftr2trk([directory,ftrfname,'_anat.mat'],refanat);

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
[wfibsmm_mni, wfibsvox_mni] = ea_map_coords(wfibsvox_anat', ...
                                            refanat, ...
                                            [directory,'inverseTransform'], ...
                                            refnorm);

wfibsmm_mni = wfibsmm_mni';
wfibsvox_mni = wfibsvox_mni';

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
if ~exist([directory,'connectomes',filesep,'dMRI'],'file')
    mkdir([directory,'connectomes',filesep,'dMRI']);
end
ea_savefibertracts([directory,'connectomes',filesep,'dMRI',filesep,ftrbase,'.mat'],[wfibsmm_mni,fibers(:,4)],idx,'mm');
ea_savefibertracts([directory,'connectomes',filesep,'dMRI',filesep,ftrbase,'_vox.mat'],[wfibsvox_mni,fibers(:,4)],idx,'vox',refnorm);

%% create normalized trackvis version
fprintf('\nExporting normalized fibers to TrackVis...\n');

[~,ftrfname]=fileparts(options.prefs.FTR_normalized);
ea_ftr2trk([directory,'connectomes',filesep,'dMRI',filesep,ftrfname]); % export normalized ftr to .trk
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
    ea_error('Please run normalization for this subject first.');
end

% determine the refimage for b0 and anat space visualization
refb0=[directory,options.prefs.b0];
refanat=[directory,options.prefs.prenii_unnormalized];

% determin the template for fiber normalization and visualization
if ismember(whichnormmethod,{'ea_normalize_spmshoot','ea_normalize_spmdartel','ea_normalize_spmnewseg'})
	refnorm=[refnorm,',2'];
end
