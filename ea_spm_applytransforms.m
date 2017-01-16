function ea_spm_applytransforms(options,usesuit)

if ~exist('usesuit','var')
    usesuit=0;
end

directory=[options.root,options.patientname,filesep];

switch options.modality
    case 1 % MR
        postops{1}=options.prefs.prenii_unnormalized;
        postops{2}=options.prefs.tranii_unnormalized;
        postops{3}=options.prefs.cornii_unnormalized;
        postops{4}=options.prefs.sagnii_unnormalized;
        gfis{1}=options.prefs.gprenii;
        gfis{2}=options.prefs.gtranii;
        gfis{3}=options.prefs.gcornii;
        gfis{4}=options.prefs.gsagnii;
        lfis{1}=options.prefs.prenii;
        lfis{2}=options.prefs.tranii;
        lfis{3}=options.prefs.cornii;
        lfis{4}=options.prefs.sagnii;
    case 2 % CT
        postops{1}=options.prefs.prenii_unnormalized;
        postops{2}=options.prefs.ctnii_coregistered;
        gfis{1}=options.prefs.gprenii;
        gfis{2}=options.prefs.gctnii;
        lfis{1}=options.prefs.prenii;
        lfis{2}=options.prefs.ctnii;
end
[postops,gfis]=ea_appendgrid(options,postops,gfis,0);

if usesuit
postops=ea_cropsuit(postops,directory);
end
% export glfiles (a bit more coarse resolution, full brain bounding box).
for pos=1:length(gfis)
    if exist([directory,postops{pos}],'file')
        nii=ea_load_nii([directory,postops{pos}]);
        
        gaussdim=abs(nii.voxsize);
        %gaussdim=abs(gaussdim(1:3)).*2;
        if mean(gaussdim>1)
            resize_img([directory,postops{pos}],gaussdim./2,nan(2,3),0);
        else
            copyfile([directory,postops{pos}],[directory,'r',postops{pos}]);
        end
        matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1}=[directory,'r',postops{pos},''];
        matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
        matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {directory};
        matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[ea_space(options),options.primarytemplate,'.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
        matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = gaussdim;
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
        clear matlabbatch jobs;
        try movefile([directory,'swr',postops{pos}],[directory,gfis{pos}]); end
    end
end

% export lfiles (fine resolution, small bounding box).
try
    for pos=1:length(lfis)
        if exist([directory,postops{pos}],'file')
            nii=ea_load_nii([directory,postops{pos}]);
            gaussdim=abs(nii.voxsize);
            %gaussdim=abs(gaussdim(1:3)).*2;
            if mean(gaussdim>1)
                resize_img([directory,postops{pos}],gaussdim./2,nan(2,3),0);
            else
                copyfile([directory,postops{pos}],[directory,'r',postops{pos}]);
            end
            
            matlabbatch{1}.spm.util.defs.comp{1}.def = {[directory,'y_ea_inv_normparams.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.push.fnames{1}=[directory,'r',postops{pos},''];
            matlabbatch{1}.spm.util.defs.out{1}.push.weight = {''};
            matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr = {directory};
            matlabbatch{1}.spm.util.defs.out{1}.push.fov.file = {[ea_space(options),'bb.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.push.preserve = 0;
            matlabbatch{1}.spm.util.defs.out{1}.push.fwhm = gaussdim;
            jobs{1}=matlabbatch;
            spm_jobman('run',jobs);
            clear matlabbatch jobs;
        end
        try movefile([directory,'swr',postops{pos}],[directory,lfis{pos}]); end
    end
end

ea_delete([directory,'rgrid.nii']);
switch options.modality
    case 1
        ea_delete([directory,'r',options.prefs.prenii_unnormalized]);
        ea_delete([directory,'r',options.prefs.tranii_unnormalized]);
        ea_delete([directory,'r',options.prefs.cornii_unnormalized]);
        ea_delete([directory,'r',options.prefs.sagnii_unnormalized]);
    case 2
        ea_delete([directory,'r',options.prefs.ctnii_coregistered]);
        ea_delete([directory,'r',options.prefs.prenii_unnormalized]);
end


function postops=ea_cropsuit(postops,directory)

[~,anatbase]=fileparts(postops{1});
% build up S for Suit:
for p=1:length(postops)
    S.resample{p}=[directory,postops{p},',1'];
end
% load affine:
load([directory,'Affine_',anatbase,'_seg1.mat']);

% load flowfield:
Vff=spm_vol([directory,'u_a_',anatbase,'_seg1.nii,1']);
% Generate grid in the space of the template
[X,Y,Z]=ndgrid(1:Vff.dim(1),1:Vff.dim(2),1:Vff.dim(3));
num_slice=Vff.dim(3);

for i=1:length(S.resample)
    [image_dir,image_name,image_ext,image_num]=spm_fileparts(S.resample{i});
    V(i)=spm_vol(S.resample{i});
    [Xm,Ym,Zm]=spmj_affine_transform(X,Y,Z,inv(Affine*V(i).mat)*Vff.mat);
    for z=1:num_slice
        Data(:,:,z)=spm_sample_vol(V(i),Xm(:,:,z),Ym(:,:,z),Zm(:,:,z),3);
    end;
    %    Data=Data.*MaskData; % Mask the image
    O(i)=V(i);
    O(i).dim=Vff.dim;
    O(i).mat=Vff.mat;
    J.images{i}{1}=fullfile(image_dir,['atemp_' image_name image_ext]);
    if (~isfield(S,'outname'))
        finalname{i}=fullfile(image_dir,['t', image_name image_ext]);
    else
        finalname{i}=S.outname{i};
    end;
    O(i).fname=J.images{i}{1};
    if (O(i).dt(1)==2)
        O(i).pinfo=[1./254*max(Data(:));0;O(i).pinfo(3)];
    end;
    spm_write_vol(O(i),Data);
    
    
    movefile([directory,postops{i}],[directory,'orig_',postops{i}]);
    movefile([directory,'atemp_',postops{i}],[directory,postops{i}]);
    
end;


function resize_img(imnames, Voxdim, BB, ismask)
%  resize_img -- resample images to have specified voxel dims and BBox
% resize_img(imnames, voxdim, bb, ismask)
%
% Output images will be prefixed with 'r', and will have voxel dimensions
% equal to voxdim. Use NaNs to determine voxdims from transformation matrix
% of input image(s).
% If bb == nan(2,3), bounding box will include entire original image
% Origin will move appropriately. Use world_bb to compute bounding box from
% a different image.
%
% Pass ismask=true to re-round binary mask values (avoid
% growing/shrinking masks due to linear interp)
%
% See also voxdim, world_bb

% Based on John Ashburner's reorient.m
% http://www.sph.umich.edu/~nichols/JohnsGems.html#Gem7
% http://www.sph.umich.edu/~nichols/JohnsGems5.html#Gem2
% Adapted by Ged Ridgway -- email bugs to drc.spm@gmail.com

% This version doesn't check spm_flip_analyze_images -- the handedness of
% the output image and matrix should match those of the input.

% Check spm version:
if exist('spm_select','file') % should be true for spm5
    spm5 = 1;
elseif exist('spm_get','file') % should be true for spm2
    spm5 = 0;
else
    error('Can''t find spm_get or spm_select; please add SPM to path')
end

spm_defaults;

% prompt for missing arguments
if ( ~exist('imnames','var') || isempty(char(imnames)) )
    if spm5
        imnames = spm_select(inf, 'image', 'Choose images to resize');
    else
        imnames = spm_get(inf, 'img', 'Choose images to resize');
    end
end
% check if inter fig already open, don't close later if so...
Fint = spm_figure('FindWin', 'Interactive'); Fnew = [];
if ( ~exist('Voxdim', 'var') || isempty(Voxdim) )
    Fnew = spm_figure('GetWin', 'Interactive');
    Voxdim = spm_input('Vox Dims (NaN for "as input")? ',...
        '+1', 'e', '[nan nan nan]', 3);
end
if ( ~exist('BB', 'var') || isempty(BB) )
    Fnew = spm_figure('GetWin', 'Interactive');
    BB = spm_input('Bound Box (NaN => original)? ',...
        '+1', 'e', '[nan nan nan; nan nan nan]', [2 3]);
end
if ~exist('ismask', 'var')
    ismask = false;
end
if isempty(ismask)
    ismask = false;
end

% reslice images one-by-one
vols = spm_vol(imnames);
for V=vols'
    % (copy to allow defaulting of NaNs differently for each volume)
    voxdim = Voxdim;
    bb = BB;
    % default voxdim to current volume's voxdim, (from mat parameters)
    if any(isnan(voxdim))
        vprm = spm_imatrix(V.mat);
        vvoxdim = vprm(7:9);
        voxdim(isnan(voxdim)) = vvoxdim(isnan(voxdim));
    end
    voxdim = voxdim(:)';

    mn = bb(1,:);
    mx = bb(2,:);
    % default BB to current volume's
    if any(isnan(bb(:)))
        vbb = world_bb(V);
        vmn = vbb(1,:);
        vmx = vbb(2,:);
        mn(isnan(mn)) = vmn(isnan(mn));
        mx(isnan(mx)) = vmx(isnan(mx));
    end

    % voxel [1 1 1] of output should map to BB mn
    % (the combination of matrices below first maps [1 1 1] to [0 0 0])
    mat = spm_matrix([mn 0 0 0 voxdim])*spm_matrix([-1 -1 -1]);
    % voxel-coords of BB mx gives number of voxels required
    % (round up if more than a tenth of a voxel over)
    imgdim = ceil(mat \ [mx 1]' - 0.1)';

    % output image
    VO            = V;
    [pth,nam,ext] = fileparts(V.fname);
    VO.fname      = fullfile(pth,['r' nam ext]);
    VO.dim(1:3)   = imgdim(1:3);
    VO.mat        = mat;
    VO = spm_create_vol(VO);
    spm_progress_bar('Init',imgdim(3),'reslicing...','planes completed');
    for i = 1:imgdim(3)
        M = inv(spm_matrix([0 0 -i])*inv(VO.mat)*V.mat);
        img = spm_slice_vol(V, M, imgdim(1:2), 1); % (linear interp)
        if ismask
            img = round(img);
        end
        spm_write_plane(VO, img, i);
        spm_progress_bar('Set', i)
    end
    spm_progress_bar('Clear');
end
% call spm_close_vol if spm2
if ~spm5
    spm_close_vol(VO);
end
if (isempty(Fint) && ~isempty(Fnew))
    % interactive figure was opened by this script, so close it again.
    close(Fnew);
end
disp('Done.')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bb = world_bb(V)
%  world-bb -- get bounding box in world (mm) coordinates

d = V.dim(1:3);
% corners in voxel-space
c = [ 1    1    1    1
    1    1    d(3) 1
    1    d(2) 1    1
    1    d(2) d(3) 1
    d(1) 1    1    1
    d(1) 1    d(3) 1
    d(1) d(2) 1    1
    d(1) d(2) d(3) 1 ]';
% corners in world-space
tc = V.mat(1:3,1:4)*c;

% bounding box (world) min and max
mn = min(tc,[],2)';
mx = max(tc,[],2)';
bb = [mn; mx];
