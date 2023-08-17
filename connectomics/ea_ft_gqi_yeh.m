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

redo=ea_prepare_dti(options);

vizz=0;

% build white matter mask
if ~exist([directory,'ttrackingmask.nii'],'file') || redo || ...
        (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    ea_gentrackingmask_brainmask(options,1)
end

basedir = [options.earoot, 'ext_libs',filesep,'dsi_studio',filesep];
DSISTUDIO = ea_getExec([basedir, 'dsi_studio'], escapePath = 1);


if options.lc.struc.ft.upsample.how==1 % internal upsampling used
    ea_roi2txt([directory,'ttrackingmask.nii'],[directory,'ttrackingmask.txt'],ea_resolve_usfactor(options.lc.struc.ft.upsample))
    maskext='.txt';
else
    maskext='.nii';
end

% build .fib.gz file
[~,ftrbase]=fileparts(options.prefs.FTR_unnormalized);
if ~exist([directory,ftrbase,'.fib.gz'],'file') || redo || ...
        (isfield(options, 'overwriteapproved') && options.overwriteapproved)
    disp('Estimating ODF / preparing GQI...');

    bvals = load([directory,options.prefs.bval]);
    if size(bvals,1)>size(bvals,2)
        bvals=bvals';
    end

    bvecs = load([directory,options.prefs.bvec]);
    if size(bvecs,1)>size(bvecs,2)
        bvecs=bvecs';
    end

    btable=[bvals;bvecs];
    ea_prepare_fib_gqi(DSISTUDIO,btable,1.2,options);

    disp('Done.');
else
    disp('.fib.gz file found, no need to rebuild.');
end

trkcmd=[DSISTUDIO,' --action=trk',...
    ' --method=0',...
    ' --source=',ea_path_helper([directory,ftrbase,'.fib.gz']),...
    ' --seed=',ea_path_helper([directory,'ttrackingmask',maskext]),...
    ' --fiber_count=', num2str(options.lc.struc.ft.dsistudio.fiber_count),...
    ' --output=',ea_path_helper([directory,ftrbase,'.mat']),...
    ' --dt_threshold=0.2',...
    ' --initial_dir=0',...
    ' --interpolation=0',...
    ' --max_length=300.0',...
    ' --min_length=10.0',...
    ' --random_seed=0',...
    ' --seed_plan=0',...
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
fibinfo = load([directory,ftrbase,'.mat']);
fibers = fibinfo.tracts';
idx = double(fibinfo.length)';
fibers = [fibers, repelem(1:numel(idx), idx)'];
clear fibinfo
b0 = spm_vol([directory,options.prefs.b0]);

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

if options.lc.struc.ft.upsample.factor == 1 % No upsampling
    % Change ZERO-BASED indexing to ONE-BASED indexing.
    fibers(:,1:3) = fibers(:,1:3) + 1;
else
    fibers = ea_resolve_usfibers(options,fibers);
end

if vizz
    figure
    thresh=700; % set to a good grey value.
    plot3(fibers(:,1),fibers(:,2),fibers(:,3),'r.')
    hold on
    b0=ea_load_nii([directory,options.prefs.b0]);
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img(:)>thresh));
    plot3(xx,yy,zz,'g.')
end

% Add index column
fibers(:,4) = repelem(1:length(idx), idx)';

ftr.ea_fibformat = '1.0';
ftr.fourindex = 1;
ftr.fibers = fibers;
ftr.idx = idx;
ftr.voxmm = 'vox';
ftr.mat = b0.mat;

disp('Saving fibers...');
save([directory,ftrbase,'.mat'],'-struct','ftr','-v7.3');
disp('Done.');

fprintf('\nGenerating trk in b0 space...\n');
ea_ftr2trk([directory,ftrbase,'.mat'], [directory,options.prefs.b0])


function ea_prepare_fib_gqi(DSISTUDIO,btable,mean_diffusion_distance_ratio,options)
directory=[options.root,options.patientname,filesep];
[~,ftrbase]=fileparts(options.prefs.FTR_unnormalized);

% source images
ea_delete([directory,'dti.src.gz']);
cmd=[DSISTUDIO,' --action=src --source=',ea_path_helper([directory,options.prefs.dti]),...
    ' --bval=',ea_path_helper([directory,options.prefs.bval])...
    ' --bvec=',ea_path_helper([directory,options.prefs.bvec])...
    ' --sort_b_table=0',...
    ' --output=',ea_path_helper([directory,'dti.src.gz'])];

if options.lc.struc.ft.upsample.how==1 % internal upsampling used
    cmd=[cmd,...
        ' --up_sampling=',num2str(factor2dsistudiofactor(ea_resolve_usfactor(options.lc.struc.ft.upsample)))];

    %% add methods dump:
    cits={
        'Dyrby, T. B., Lundell, H., Burke, M. W., Reislev, N. L., Paulson, O. B., Ptito, M., & Siebner, H. R. (2013). Interpolation of diffusion weighted imaging datasets. NeuroImage, 103(C), 1?12. http://doi.org/10.1016/j.neuroimage.2014.09.005'
        'Yeh, F.-C., Wedeen, V. J., & Tseng, W.-Y. I. (2010). Generalized q-Sampling Imaging. IEEE Transactions on Medical Imaging, 29(9), 1626?1635. http://doi.org/10.1109/TMI.2010.2045126'
        };
    ea_methods(options,['Raw diffusion data was upsampled using bspline-interpolation with a factor of ',num2str(factor2dsistudiofactor(ea_resolve_usfactor(options.lc.struc.ft.upsample))),' following the concept described in Dyrby et al. 2014 as implemented in DSI Studio (http://dsi-studio.labsolver.org/; Yeh et al. 2010).'],cits);

    maskext='.txt';
else
    maskext='.nii';
end

err=ea_runcmd(cmd);

if err || ~exist([directory,'dti.src.gz'],'file')
    ea_warning('DSI studio failed to generate .src file. Using Matlab code instead.');
    ea_create_dsistudio_src([directory,options.prefs.dti],[directory,'dti.src']);
end

% create .fib file
cmd=[DSISTUDIO,' --action=rec --source=',ea_path_helper([directory,'dti.src.gz']),...
    ' --mask=',ea_path_helper([directory,'ttrackingmask',maskext])...
    ' --method=4',...
    ' --param0=1.25',...
    ' --num_fiber=5',...
    ' --odf_order=8'];

err=ea_runcmd(cmd);
ea_delete([directory,'dti.src.gz']);

if err
    disp('Reconstruction from command line failed. Reattempting inside Matlab.');
    % do it the matlab way
    res=ea_gqi_reco([directory,options.prefs.dti],btable,mean_diffusion_distance_ratio,options);
    delete([directory,ftrbase,'.fib']);
    save([directory,ftrbase,'.fib'],'-struct','res','-v4');
    if ~exist([directory,ftrbase,'.fib'],'file')
       ea_error(['It seems like the .fib file created for ',directory,' is too large for DSI studio to handle. We are working on a solution..']);
    end
    gzip([directory,ftrbase,'.fib']);
    ea_delete([directory,ftrbase,'.fib']);
else
    di=dir([directory,'dti.src.gz*.fib.gz']);
    if length(di)>1
        ea_error('Too many .fib.gz files present in folder. Please delete older files');
    end
    movefile([directory,di(1).name],[directory,ftrbase,'.fib.gz']);
end


%% add methods dump:
cits={
    'Yeh, F.-C., Wedeen, V. J., & Tseng, W.-Y. I. (2010). Generalized q-sampling imaging. IEEE Transactions on Medical Imaging, 29(9), 1626?1635. http://doi.org/10.1109/TMI.2010.2045126'
    'Ashburner, J., & Friston, K. J. (2005). Unified segmentation., 26(3), 839?851. http://doi.org/10.1016/j.neuroimage.2005.02.018'
    };
ea_methods(options,['A whole-brain fiber-set was estimated based using the Generalized q-sampling imaging (GQI) approach (Yeh 2010) as implemented in DSI-Studio (http://dsi-studio.labsolver.org).',...
    ' GQI is a model-free method that calculates the orientational distribution of the density of diffusing water.',...
    ' Fibers were sampled within a white-matter mask that was estimated using the anatomical acquisition by applying the Unified Segmentation approach (Ashburner 2005) as implemented in ',spm('ver'),'. This mask was linearly co-registered to the b0-weighted series.'],cits);



function res=ea_gqi_reco(filename,b_table,mean_diffusion_distance_ratio,options)
% Direct GQI reconstruction from huge image data
% You may need to include find_peak.m to run these codes.
% parameters:
%  filename: the filename of the image volume
%  file_type: the pixel format of the image, can be 'int8', 'int16', 'int32', 'single', or 'double'
%  pixel_size: the size of the pixel in bytes.
%  b_table: the b-table matrix with size of 4-by-d, where d is the number of the diffusion weighted images
%          b(1,:) stores the b-value, whereas b(2:4,:) stores the grandient vector
%  mean_diffusion_distance_ratio: check out GQI reconstruction for detail. Recommended value=1.2
%
% example:
%  gqi_reco('2dseq','int32',4,b_table,1.0);

load([options.earoot,'connectomics',filesep,'ea_gqi_odf8.mat']);

nii=ea_load_nii(filename);
dim = size(nii(1).img);
dif = length(nii);

for d=1:3
    ptm=[1,1,1,1];
    m1=nii(1).mat*ptm';
    ptm(d)=2;
    m2=nii(1).mat*ptm';
    res.voxel_size(d) = ea_pdist([m1,m2]');
end

fa0 = zeros(dim);
fa1 = zeros(dim);
fa2 = zeros(dim);
index0 = zeros(dim);
index1 = zeros(dim);
index2 = zeros(dim);

reco_temp = zeros(dim(1),dim(2),dif);
plane_size = dim(1)*dim(2);

% GQI reconstruction matrix A
l_values = sqrt(b_table(1,:)*0.01506);
b_vector = b_table(2:4,:).*repmat(l_values,3,1);
A = sinc(odf_vertices'*b_vector*mean_diffusion_distance_ratio/pi);

%f =fopen(filename);
max_dif = 0;

for z = 1:dim(3)
    for d = 1:dif
        %fseek(f,((z-1)*plane_size+(d-1)*dim(1)*dim(2)*dim(3))*pixel_size,'bof');
        reco_temp(:,:,d) = nii(d).img(:,:,z);
    end
    for x = 1:dim(1)
        for y = 1:dim(2)
            ODF=A*reshape(reco_temp(x,y,:),[],1);
            ODF=ODF(:);
            p = ea_find_peak(ODF,odf_faces);
            max_dif = max(max_dif,mean(ODF));
            min_odf = min(ODF(:));
            fa0(x,y,z) = ODF(p(1))-min_odf;
            index0(x,y,z) = p(1)-1;
            if length(p) > 1
                fa1(x,y,z) = ODF(p(2))-min_odf;
                index1(x,y,z) = p(2)-1;
            end
            if length(p) > 2
                fa2(x,y,z) = ODF(p(3))-min_odf;
                index2(x,y,z) = p(3)-1;
            end
        end
    end
end
res.fa0 = fa0/max_dif;
res.fa1 = fa1/max_dif;
res.fa2 = fa2/max_dif;
res.fa0 = reshape(res.fa0,1,[]);
res.fa1 = reshape(res.fa1,1,[]);
res.fa2 = reshape(res.fa2,1,[]);
res.index0 = reshape(index0,1,[]);
res.index1 = reshape(index1,1,[]);
res.index2 = reshape(index2,1,[]);
res.dimension = dim;
res.odf_faces=odf_faces;
res.odf_vertices=odf_vertices;

function ofactor=factor2dsistudiofactor(factor)
switch factor
    case 2
        ofactor=1; % DSI studio codes upsampling factor of 2 with 1
    case 4
        ofactor=3; % DSI studio codes upsampling factor of 4 with 3
end

function p = ea_find_peak(odf,odf_faces)
is_peak = odf;
odf_faces = odf_faces + 1;
odf_faces = odf_faces - (odf_faces > length(odf))*length(odf);
is_peak(odf_faces(1,odf(odf_faces(2,:)) >= odf(odf_faces(1,:)) | ...
    odf(odf_faces(3,:)) >= odf(odf_faces(1,:)))) = 0;
is_peak(odf_faces(2,odf(odf_faces(1,:)) >= odf(odf_faces(2,:)) | ...
    odf(odf_faces(3,:)) >= odf(odf_faces(2,:)))) = 0;
is_peak(odf_faces(3,odf(odf_faces(2,:)) >= odf(odf_faces(3,:)) | ...
    odf(odf_faces(1,:)) >= odf(odf_faces(3,:)))) = 0;
[values,ordering] = sort(-is_peak);
p = ordering(values <= 0);
