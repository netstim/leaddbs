function ea_normalize_fibers(options)
% uses map_coords function by Ged Ridgway (see below)
directory=[options.root,options.patientname,filesep];

% create (unnormalized) trackvis version
disp('Exporting to TrackVis');
[~,ftrfname]=fileparts(options.prefs.FTR_unnormalized);
try
%     if ~exist([directory,ftrfname,'.trk'],'file')
        reftemplate=[directory,options.prefs.b0];
        dnii=ea_load_nii(reftemplate);
        niisize=size(dnii.img); % get dimensions of reference template.

        specs.origin=[0,0,0];
        specs.dim=niisize;
        specs.vox=ftr.vox;
        specs.affine=dnii.mat;

        ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
%     end
end
disp('Done.');

if ~exist([directory,'y_ea_inv_normparams.nii'],'file');
    ea_error('Please run a compatible normalization of the preoperative MRI-volume first. Final (inverse) normalization parameters should be stored as y_ea_inv_normparams.nii inside of the subject folder.');
end

vizz=0; % turn this value to 1 to visualize fiber normalization (option for debugging only, this will drastically slow down the process).
cleanse_fibers=0; % deletes everything outside the white matter of the template.

% check which normalization routine has been used..
% if dartel was used, we need to coregister c2 of b0 and rc2 of anat (since
% deformation fields were estimated for the rc* files and not the native
% anat file.
[options.prefs.b0,options.prefs.prenii_unnormalized]=ea_checkdartelused(options);

%% normalize fibers

% get affinematrix from b0 to preop mri
Vfirst=spm_vol([directory,options.prefs.b0,',1']);
Vsecond=spm_vol([directory,options.prefs.prenii_unnormalized,',1']);
x=spm_coreg(Vfirst,Vsecond);
affinematrix1=Vsecond.mat\spm_matrix(x(:)')*Vfirst.mat;

b0=ea_load_nii([directory,options.prefs.b0]);

ysize=size(b0.img,2)+1;

ftr = load([directory,options.prefs.FTR_unnormalized]);

if ~exist([options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6_hires.nii'],'file')
    gunzip([options.earoot,'templates',filesep,'dartel',filesep,'create_mni_darteltemplate',filesep,'dartelmni_6_hires.nii.gz'],...
        [options.earoot,'templates',filesep,'dartel']);
end

reftemplate=[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6_hires.nii,2'];
Vmni=spm_vol(reftemplate);
mnimask=spm_read_vols(Vmni);
mnimask=mnimask>0.01;

if vizz
    figure('color','w');
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
    title('Anat space');
    [xx,yy,zz]=ind2sub(size(anat.img),find(anat.img>max(anat.img(:))/3));
    plot3(xx(1:1000:end),yy(1:1000:end),zz(1:1000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
    % plot MNI
    mni=ea_load_nii([options.earoot,'templates',filesep,'mni_hires.nii']);
    subplot(1,3,3);
    title('MNI space');
    [xx,yy,zz]=ind2sub(size(mni.img),find(mni.img>max(mni.img(:))/3));
    % transpose to mm
    XYZ=[xx,yy,zz,ones(length(xx),1)]';
    XYZ=mni.mat*XYZ;
    plot3(XYZ(1,1:10000:end),XYZ(2,1:10000:end),XYZ(3,1:10000:end),'.','color',[0.9598    0.9218    0.0948]);
    axis vis3d off tight equal;
    hold on
end

ea_dispercent(0,'Normalizing fibers');
numfibs=length(ftr.curveSegCell);

ynii=nifti([directory,'y_ea_inv_normparams.nii']);
P = [repmat([directory,'y_ea_inv_normparams.nii'],3,1),[',1,1';',1,2';',1,3']];
Vnii = spm_vol(P);
wfibs=cell(length(ftr.curveSegCell),1);
deletefibers=[];
for fib=1:numfibs

    ea_dispercent(fib/numfibs);

    %% transpose from freiburg to spm notation.
    %wfibs{fib}=[ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,2),ftr.curveSegCell{fib}(:,3),ones(length(ftr.curveSegCell{fib}),1)];
    wfibs{fib}=[ftr.curveSegCell{fib}(:,2),ysize-ftr.curveSegCell{fib}(:,1),ftr.curveSegCell{fib}(:,3),ones(length(ftr.curveSegCell{fib}),1)];
    if vizz
       thisfib=wfibs{fib}';
        subplot(1,3,1)
        plot3(thisfib(1,:),thisfib(2,:),thisfib(3,:),'-','color',[0.1707    0.2919    0.7792]);
    end

    %% first apply affine transform from b0 to prenii
    wfibs{fib}=affinematrix1*wfibs{fib}';
    if vizz
        thisfib=wfibs{fib};
        subplot(1,3,2)
        plot3(thisfib(1,:),thisfib(2,:),thisfib(3,:),'-','color',[0.1707    0.2919    0.7792]);
    end
    %% -> coordinates are now in voxel-space of single subject anat file.

    %% map from prenii voxelspace to mni millimeter space

    wfibs{fib} = vox2mm_mni(wfibs{fib},Vnii,ynii)';
    if vizz
        thisfib=wfibs{fib}';
        subplot(1,3,3)
        plot3(thisfib(1,:),thisfib(2,:),thisfib(3,:),'-','color',[0.1707    0.2919    0.7792]);
    end

    %% map from mni millimeter space to mni voxel space (only needed for trackvis convertion and cleansing fibers).
    wfibsvox{fib}=[wfibs{fib},ones(size(wfibs{fib},1),1)]';
    wfibsvox{fib}=Vmni(1).mat\wfibsvox{fib};
    wfibsvox{fib}=wfibsvox{fib}(1:3,:)';
    wfibsvox{fib}=[wfibsvox{fib}(:,1),wfibsvox{fib}(:,2),wfibsvox{fib}(:,3)];

    %% cleansing fibers..
    if cleanse_fibers % delete anything too far from wm.

        todelete=~mnimask(sub2ind(size(mnimask),round(wfibsvox{fib}(:,1)),round(wfibsvox{fib}(:,2)),round(wfibsvox{fib}(:,3))));

        if all(todelete) % all fibers outside WM
            deletefibers=[deletefibers,fib];
        else
            wfibs{fib}(todelete,:)=[];
            wfibsvox{fib}(todelete,:)=[];
        end

    end

    %% cleanup
    %wfibs{fib}=wfibs{fib}(:,1:3);
   if vizz; drawnow; end
end

wfibs(deletefibers)=[]; % delete fibers that were in total outside WM
wfibsvox(deletefibers)=[]; % delete fibers that were in total outside WM

ea_dispercent(100,'end');

wfibsvox=wfibsvox';
nftr.normalized_fibers_mm=wfibs; clear wfibs
nftr.normalized_fibers_vox=wfibsvox; clear wfibsvox
if isfield(ftr,'curveD')
    nftr.curveD=ftr.curveD;
end
nftr.trackParam=ftr.trackParam;
nftr.user=ftr.user;
nftr.vox=Vmni.mat(logical(eye(4))); nftr.vox=nftr.vox(1:3)';
disp('Saving files...');
save([directory,options.prefs.FTR_normalized],'-struct','nftr','-v7.3');
%save([directory,'vox_',options.prefs.FTR_normalized],'normalized_fibers_vox');
disp('Done.');

%% create trackvis version
disp('Creating TrackVis version...');
try
    reftemplate=[options.earoot,'templates',filesep,'dartel',filesep,'dartelmni_6_hires.nii'];
    dnii=ea_load_nii(reftemplate);
    niisize=size(dnii(1).img); % get dimensions of reference template.

    specs.origin=[0,0,0];
    specs.dim=niisize;
    specs.vox=ftr.vox;
    specs.affine=dnii(1).mat;

    [~,ftrfname]=fileparts(options.prefs.FTR_normalized);
    ea_ftr2trk(ftrfname,directory,specs,options); % export normalized ftr to .trk
end
% delete([directory,'vox_',options.prefs.FTR_normalized]);

disp('Done.');


function [useb0,useanat]=ea_checkdartelused(options)
directory=[options.root,options.patientname,filesep];

[whichnormmethod,tempfile]=ea_whichnormmethod(directory);
switch whichnormmethod
    case 'ea_normalize_ants'
        ea_error('ANTs normalization is not supported for Fibers normalization right now.');
    case 'ea_normalize_spmdartel'
        dartelused=1;
    otherwise
        dartelused=0;
end

if dartelused
    % segment b0.
    if ~exist([directory,'c2',options.prefs.b0],'file');
        disp('Segmenting B0 file for DARTEL import space coregistration...');
        ea_newseg(directory,options.prefs.b0,0,options);
        delete([directory,'c4',options.prefs.b0]);
        delete([directory,'c5',options.prefs.b0]);
        disp('Done.');
    end

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
        cfg_util('run',{matlabbatch});
        clear matlabbatch
    end
end

if ~dartelused
    useb0=options.prefs.b0;
    useanat=options.prefs.prenii_unnormalized;
else
    useb0=['c2',options.prefs.b0];
    useanat=['rc2',options.prefs.prenii_unnormalized];
end


function coord = vdox2mm_mni(coord,Vnii,ynii)
ixs = double(coord(1:3, :));

% old method
for i = 1:3
coord(i,:)=ynii.dat(sub2ind(size(ynii.dat),ixs(1,:)',ixs(2,:)',ixs(3,:)',ones(size(ixs,2),1),repmat(i,size(ixs,2),1)));
end

function coord = vox2mm_mni(coord, Vnii,ynii)
% new method
ixs_new = double(coord(1:3, :));
coord=[spm_sample_vol(Vnii(1),ixs_new(1,:),ixs_new(2,:),ixs_new(3,:),1);
    spm_sample_vol(Vnii(2),ixs_new(1,:),ixs_new(2,:),ixs_new(3,:),1);
    spm_sample_vol(Vnii(3),ixs_new(1,:),ixs_new(2,:),ixs_new(3,:),1)];


% %old method
% ixs_old = round(coord(1:3, :));
% coord=zeros(3,size(coord,2));
% for i = 1:3
% coord(i,:)=ynii.dat(sub2ind(size(ynii.dat),ixs_old(1,:)',ixs_old(2,:)',ixs_old(3,:)',ones(size(ixs_old,2),1),repmat(i,size(ixs_old,2),1)));
% end

