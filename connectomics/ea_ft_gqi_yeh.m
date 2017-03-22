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
ea_prepare_dti(options)

vizz=0;

load([options.root,options.patientname,filesep,options.prefs.bval]);
[~,bvfname]=fileparts(options.prefs.bval);

bvals=eval(bvfname);
if size(bvals,1)>size(bvals,2)
    bvals=bvals';
end

load([options.root,options.patientname,filesep,options.prefs.bvec]);
[~,bvfname]=fileparts(options.prefs.bvec);
bvecs=eval(bvfname);
if size(bvecs,1)>size(bvecs,2)
    bvecs=bvecs';
end
btable=[bvals;bvecs];

% build white matter mask

if ~exist([options.root,options.patientname,filesep,'ttrackingmask.nii'],'file');
ea_gentrackingmask(options,1)
end

basedir = [options.earoot, 'ext_libs',filesep,'dsi_studio',filesep];
if ismac
    dsistudio = [basedir,'mac',filesep, 'dsi_studio.app',filesep,'Contents',filesep,'MacOS',filesep,'dsi_studio'];
elseif isunix
    ea_libs_helper([basedir, 'linux']);
    dsistudio = [basedir, 'linux',filesep,'dsi_studio'];
elseif ispc
    dsistudio = [basedir, 'win',filesep,'dsi_studio.exe'];
end

% build .fib.gz file
[~,ftrbase]=fileparts(options.prefs.FTR_unnormalized);
if ~exist([options.root,options.patientname,filesep,ftrbase,'.fib.gz'],'file')
    disp('Estimating ODF / preparing GQI...');
    ea_prepare_fib_gqi(dsistudio,btable,1.2,options);

    disp('Done.');
else
    disp('.fib.gz file found, no need to rebuild.');
end


cmd=[dsistudio,' --action=trk --source=',ea_path_helper([options.root,options.patientname,filesep,ftrbase,'.fib.gz']),...
    ' --method=0',...
    ' --seed=',ea_path_helper([options.root,options.patientname,filesep,'ttrackingmask.nii']),...
    ' --fiber_count=50000',...
    ' --output=',ea_path_helper([options.root,options.patientname,filesep,ftrbase,'.mat'])];

err=ea_submitcmd(cmd);
if err
    ea_error(['Fibertracking with dsi_studio failed (error code=',num2str(err),').']);
end

% now store tract in lead-dbs format
ea_dispercent(0,'Converting fibers');
fibinfo=load([options.root,options.patientname,filesep,ftrbase,'.mat']);
fibers=fibinfo.tracts;
idx=fibinfo.length';
clear fibinfo
fibers=fibers';
b0=spm_vol([options.root,options.patientname,filesep,options.prefs.b0]);

% default is a x- and y-flip as per communication with fang-cheng yeh
% 07-13-2016

if b0.mat(1)>0
    % flip x
    fibers(:,1)=(b0.dim(1))-fibers(:,1);
end
if b0.mat(6)>0
    %flip y
    fibers(:,2)=(b0.dim(2))-fibers(:,2);
end
if b0.mat(11)<0
    %flip z
    fibers(:,3)=(b0.dim(3))-fibers(:,3);
end

if vizz
    figure
    thresh=700; % set to a good grey value.
    plot3(fibers(:,1),fibers(:,2),fibers(:,3),'r.')
    hold on
    b0=ea_load_nii([options.root,options.patientname,filesep,options.prefs.b0]);
    [xx,yy,zz]=ind2sub(size(b0.img),find(b0.img(:)>thresh));
    plot3(xx,yy,zz,'g.')
end

clear length
idxv=zeros(size(fibers,1),1);
lid=1; cnt=1;
for id=idx'
    ea_dispercent(cnt/length(idx));
    idxv(lid:lid+id-1)=cnt;
    lid=lid+id;
    cnt=cnt+1;
end
ea_dispercent(1,'end');

fibers=[fibers,idxv];

ftr.fourindex=1;
ftr.ea_fibformat='1.0';
ftr.fibers=fibers;
ftr.idx=idx;
disp('Saving fibers...');
save([options.root,options.patientname,filesep,ftrbase,'.mat'],'-struct','ftr','-v7.3');
disp('Done.');


function ea_prepare_fib_gqi(dsistudio,btable,mean_diffusion_distance_ratio,options)
[~,ftrbase]=fileparts(options.prefs.FTR_unnormalized);

if exist([options.root,options.patientname,filesep,ftrbase,'.fib.gz'],'file')
   disp('.fib.gz file already present, no need to rebuild.');
   return
end

% try the DSI-studio way (faster):

% source images
cmd=[dsistudio,' --action=src --source=',ea_path_helper([options.root,options.patientname,filesep,options.prefs.dti]),...
    ' --bval=',ea_path_helper([options.root,options.patientname,filesep,options.prefs.bval])...
    ' --bvec=',ea_path_helper([options.root,options.patientname,filesep,options.prefs.bvec])...
    ' --output=',ea_path_helper([options.root,options.patientname,filesep,'dti.src.gz'])];

err=ea_submitcmd(cmd);

if err
    warning(['Sourcing from command line with dsi_studio failed (error code=',num2str(err),').']);
end

% create .fib file
cmd=[dsistudio,' --action=rec --source=',ea_path_helper([options.root,options.patientname,filesep,'dti.src.gz']),...
    ' --mask=',ea_path_helper([options.root,options.patientname,filesep,'ttrackingmask.nii'])...
    ' --method=4',...
    ' --param0=1.25'];

err=ea_submitcmd(cmd);
delete([options.root,options.patientname,filesep,'dti.src.gz']);
if err
    warning(['Reconstruction from command line with dsi_studio failed (error code=',num2str(err),').']);
end

di=dir([options.root,options.patientname,filesep,'dti.src.gz*.fib.gz']);
if length(di)>1
    ea_error('Too many .fib.gz files present in folder. Please delete older files');
end
movefile([options.root,options.patientname,filesep,di(1).name],[options.root,options.patientname,filesep,ftrbase,'.fib.gz']);

if ~exist([options.root,options.patientname,filesep,ftrbase,'.fib.gz'],'file');
    disp('Reconstruction from command line failed. Reattempting inside Matlab.');

    % do it the matlab way
    res=ea_gqi_reco([options.root,options.patientname,filesep,options.prefs.dti],btable,mean_diffusion_distance_ratio,options);
    save([options.root,options.patientname,filesep,ftrbase,'.fib'],'-struct','res','-v4');
    gzip([options.root,options.patientname,filesep,ftrbase,'.fib']);
    ea_delete([options.root,options.patientname,filesep,ftrbase,'.fib']);

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
    res.voxel_size(d) = pdist([m1,m2]');
end

fa0 = zeros(dim);
fa1 = zeros(dim);
fa2 = zeros(dim);
index0 = zeros(dim);
index1 = zeros(dim);
index2 = zeros(dim);

reco_temp = zeros(dim(1),dim(2),dif);
plane_size = dim(1)*dim(2);

% GQI reconstruciton matrix A
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
            p = ea_find_peak(ODF,odf_faces);
            max_dif = max(max_dif,mean(ODF));
            min_odf = min(ODF);
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
