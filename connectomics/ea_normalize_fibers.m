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
        specs.vox=dnii.voxsize;
        specs.affine=dnii.mat;

        ea_ftr2trk(ftrfname,directory,specs); % export normalized ftr to .trk
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
    hold on

    % plot anat, voxel space
    anat=ea_load_nii(refanat);
    subplot(1,3,2);
    title('anat space');
    [xx,yy,zz]=ind2sub(size(anat.img),find(anat.img>max(anat.img(:))/3));
    plot3(xx(1:1000:end),yy(1:1000:end),zz(1:1000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on

    % plot MNI, world space
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

display(sprintf('\nNormalizing fibers...'));

%% Normalize fibers

%% map from b0 voxel space to anat mm and voxel space
display(sprintf('\nMapping from b0 to anat...'));
ea_coreg2images(options,refb0,refanat,[options.root,options.patientname,filesep,'tmp.nii'],{},1);
[~, mov] = fileparts(options.prefs.b0);
[~, fix] = fileparts(options.prefs.prenii_unnormalized);
[~, wfibsvox_anat] = ea_map_coords(fibers(:,1:3)', ...
                                   refb0, ...
                                   [directory, mov, '2', fix, '.mat'], ...
                                   refanat, ...
                                   options.coregmr.method);
wfibsvox_anat = wfibsvox_anat';

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
display(sprintf('\nMapping from anat to mni...'));
[wfibsmm_mni, wfibsvox_mni] = ea_map_coords(wfibsvox_anat', ...
                                            refanat, ...
                                            [directory,'y_ea_inv_normparams.nii'], ...
                                            refnorm);

wfibsmm_mni = wfibsmm_mni';
wfibsvox_mni = wfibsvox_mni';

mniaffine=spm_get_space(refnorm);

display(sprintf('\nNormalization done.'));

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
    specs.vox=dnii.voxsize;
    specs.affine=dnii.mat;

    [~,ftrfname]=fileparts(options.prefs.FTR_normalized);
    ea_ftr2trk(ftrfname,[directory,'connectomes',filesep,'dMRI',filesep],specs); % export normalized ftr to .trk
    disp('Done.');
end




%% add methods dump:
cits={
    'Horn, A., Ostwald, D., Reisert, M., & Blankenburg, F. (2014). The structural-functional connectome and the default mode network of the human brain. NeuroImage, 102 Pt 1, 142?151. http://doi.org/10.1016/j.neuroimage.2013.09.069'
    'Horn, A., & Kühn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'
    'Horn, A., & Blankenburg, F. (2016). Toward a standardized structural-functional group connectome in MNI space. NeuroImage, 124(Pt A), 310?322. http://doi.org/10.1016/j.neuroimage.2015.08.048'
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
