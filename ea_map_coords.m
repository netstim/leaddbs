function [XYZ_dest_mm, XYZ_dest_vx] = ea_map_coords(varargin)
% Coordinates mapping, support SPM, ANTs and FSL linear and non-linear
% transformation.
%
% Please use this function for any point transform within the Lead Suite
% environment. It should automatically detect whether to use ANTs, SPM or
% FSL based transforms for the respective subject.
%
% Parameters:
%     XYZ_src_vx: coordinates to be transformed, size: 3*N or 4*N
%     src: the image defining the space of the src coordinates
%     transform: transformation from src image to dest image, optional
%     dest: the image defining the space of the dest coordinates, optional
%     transformmethod: the transformation method used, case insensitive, optional
%
%     useinverse: use the inverse of the transformation or not, only needed
%                 by ANTs in the manual mode (calling outside LEAD with
%                 explicitly specified transformation file), optional
%
% Output:
%     XYZ_dest_mm: mm coordinates in dest image. size: 3*N
%     XYZ_dest_vx: vox coordinates in dest image. size: 3*N
%
% To map the 'vox' coords in src to the 'mm' coords in src:
%     XYZ_mm = ea_map_coords(XYZ_vx, src);
%
% To map the 'mm' coords in src to the 'vox' coords in src:
%     [~, XYZ_vox] = ea_map_coords(XYZ_mm, src);
%
% To map the 'vox' coords in src to the 'mm' coords in dest:
%     XYZ_dest_mm = ea_map_coords(XYZ_src_vx, src, transform);
%
% To map the 'vox' coords in src to the 'mm' coords in dest, LeadDBS environment:
%     XYZ_dest_mm = ea_map_coords(XYZ_vx, src, 'SUBJECT_PATH/inverseTransform', dest);
%
% To map the 'vox' coords in src to the 'mm' coords in dest with transformmethod applied:
%     XYZ_dest_mm = ea_map_coords(XYZ_vx, src, transform, dest, transformmethod);
%
% 'transform' is the tranformation generated by registrating src to dest
% (can also use the inverse transformation in some cases like ANTs or FLIRT)
% It can be:
%
%     4*4 affine matrix (vox to mm affine)
%
%     4*4 affine matrix as mat file (vox to mm affine)
%     transformmethod: 'AFFINE'
%
%     1*N transformation from spm_coreg(dest, src)
%
%     1*N transformation as mat file from spm_coreg(dest, src)
%     transformmethod: 'COREG'
%
%     '*.mat' file from SPM coreg.estwrite (vox to mm affine), check ea_spm_coreg
%     transformmethod: 'SPM'
%
%     '*.mat' file from ANTs (Linear) (mm to mm affine), check ea_ants
%     transformmethod: 'ANTS'
%
%     '*.mat' file from FSL (FLIRT) (mm to mm afine), check ea_flirt
%     transformmethod: 'FSL' or 'FLIRT'
%
%     old '*_sn.mat' file from SPM (DCT structure)
%
%     'subjDir/forwardTransform' or 'subjDir/inverseTransform' stub.
%     In such case, real transformation file for ANTs, FSL or SPM will be
%     automatically detected (in the subject's folder). Thus you can
%     use an unified call to do the mapping, no need to check the
%     normalization method beforehand. Check the comments below for more
%     informations.
%
%     'y_*.nii' or 'iy_*.nii' file from SPM
%
%     '*.nii', '*.nii.gz' or '*.h5' file from ANTs (Non-linear)
%     transformmethod: 'ANTS'
%
%     '*.nii' or '*.nii.gz' files from FNIRT
%     transformmethod: 'FSL' or 'FNIRT'
%
% For SPM and ANTs, to map the coords in src image to the coords in dest
% image (the registration was done using src image as moving image and dest
% image as fixed image), the INVERSE version of the transformation should
% be used.
%
% For FSL, to map the coords in src image to the coords in dest image
% (still, the registration was done using src image as moving image and
% dest image as fixed image), the direct warp field is used as in the
% official document. But this way has severe performance issue since it
% internally inverts the warp field for each point in each iteration. To
% solve this problem, here we make a modified version of 'img2imgcoord',
% which can also use the inverse version of the warp filed. Thus the coords
% mapping is extremely speeded up. So it is recommended here to use the
% inverse version of the warp field if you want to do the coords mapping
% manually.
%
% If the registration was done by ANTs, and you want to manually do the
% coords mapping (outside the LEAD environment):
%    XYZ_dest_mm = ea_map_coords(XYZ_src_vx, src, transform, dest, 'ANTS');
% The 'tranform' here should be the deformation field file generated by
% registering src to dest (not the inverse one), 'useinverse' is set to 1
% by default. Set 'useinverse' only if you really know what you are doing.


if nargin < 2
    error('Must specify at least coords and src!')
end

XYZ_src_vx=varargin{1};

src=varargin{2};

if nargin >= 3
    transform=varargin{3};
end

if nargin >= 4
    dest=varargin{4};
end

% Check input coordinates, XYZ_vx should be column vector: 3*N or 4*N
if size(XYZ_src_vx, 1) == 3
    % make homogeneous
    XYZ_src_vx = [XYZ_src_vx; ones(1,size(XYZ_src_vx, 2))];
elseif size(XYZ_src_vx, 1) ~= 4
    error('Coord array must have 3 or 4 rows: [x;y;z] or [x;y;z;1]')
end

% srcvx to/from srcmm only
if nargin == 2
    if nargout == 1
        XYZ_dest_mm = ea_get_affine(src) * XYZ_src_vx;
    elseif nargout == 2
        % if input coords are actually in world space, then output the
        % voxel space
        XYZ_dest_mm = varargin{1};
        dest = src;
        XYZ_dest_vx = ea_get_affine(dest) \ XYZ_dest_mm;
        XYZ_dest_vx = XYZ_dest_vx(1:3,:);
    end
    XYZ_dest_mm = XYZ_dest_mm(1:3,:);
    transform = []; % finish mapping, set to empty
end

% transformation specified
if ~isempty(transform)

    % transformation is a variable, LINEAR case
    if ~ischar(transform)

        if isequal(size(transform),[4 4]) % 4*4 affine matrix supplied
            XYZ_dest_mm = transform * XYZ_src_vx;

        elseif size(transform,1) == 1 % 1*N return value from spm_coreg(dest, src) suppplied
            XYZ_dest_mm = spm_matrix(transform(:)')\ea_get_affine(src)*XYZ_src_vx;

        else
            error('Improper or unsuported transform specified!');
        end

    % DCT structure from old SPM code, NON-LINEAR case
    elseif ~isempty(regexp(transform, 'sn\.mat$', 'once'))

        % DCT sn structure
        XYZ_dest_mm = srcvx2destmm_sn(XYZ_src_vx, transform);

    % mat file supplied, LINEAR case
    elseif ~isempty(regexp(transform, '\.mat$', 'once'))

        % Need to differentiate  the transformation type
        if nargin >= 5
            normMethod = upper(varargin{5});
        else
            normMethod = 'FALLBACK';
        end


        if strcmp(normMethod, 'AFFINE') % file is  4*4 affine matrix
                transform = load(transform);
                varname = fieldnames(transform);
                transform = transform.(varname{1});
                XYZ_dest_mm = transform * XYZ_src_vx;

        elseif strcmp(normMethod, 'COREG') % file is 1*N spm_coreg return value
                transform = load(transform);
                varname = fieldnames(transform);
                transform = transform.(varname{1});
                XYZ_dest_mm = spm_matrix(transform(:)')\ea_get_affine(src)*XYZ_src_vx;

        elseif contains(normMethod, 'SPM') && ~contains(normMethod, 'HYBRID') % Registration done by SPM (ea_spm_coreg)
                % fuzzy match, transform can be specified as *_spm.mat or simply *.mat
                if ~isfile(transform) % *.mat not specified
                    transform = regexprep(transform, '\.mat$', '-spm.mat');
                    if ~isfile(transform)
                        error('Transformation file not detected!');
                    end
                end

                transform = load(transform, 'spmaffine');
                transform = transform.spmaffine;
                XYZ_dest_mm = transform * XYZ_src_vx;

         elseif contains(normMethod, 'ANTS') % Registration done by ANTs (ea_ants)
                % fuzzy match, transform can be specified as *_ants.mat or simply *.mat
                if isfile(transform) % *.mat present
                    useinverse = 1;
                else % *.mat not present
                    transform = regexprep(transform, '\.mat$', '-ants.mat');
                    if isfile(transform) % *_ants.mat present
                        useinverse = 1; % Registration was done from src to dest, so we need to use inverse here
                    else
                        % Check if inverse transformation exist
                        formSpace = regexp(transform, '(?<=_from-)[a-zA-Z]+', 'match', 'once');
                        toSpace = regexp(transform, '(?<=_to-)[a-zA-Z]+', 'match', 'once');
                        transform = strrep(transform, ['_from-', formSpace], ['_from-', toSpace]);
                        transform = strrep(transform, ['_to-', toSpace], ['_to-', formSpace]);
                        useinverse = 0;

                        if ~isfile(transform) % Inverse transformation not present
                            error('Transformation file not detected!');
                        end
                    end
                end

                % vox to mm, ANTs takes mm coords as input
                XYZ_src_mm = ea_get_affine(src)*XYZ_src_vx;

                % RAS to LPS, ANTs (ITK) use LPS coords
                XYZ_src_mm(1,:)=-XYZ_src_mm(1,:);
                XYZ_src_mm(2,:)=-XYZ_src_mm(2,:);

                % apply transform, need transpose becuase ANTs prefer N*3 like row vector
                try
                    XYZ_dest_mm = ea_antspy_apply_transforms_to_points(XYZ_src_mm(1:3,:)', transform, useinverse)';
                catch
                    ea_cprintf('CmdWinWarnings', 'Failed to run ANTsPy! Fallback to ANTs...\n');
                    XYZ_dest_mm = ea_ants_apply_transforms_to_points(XYZ_src_mm(1:3,:)', transform, useinverse)';
                end

                % LPS to RAS, restore to RAS coords
                XYZ_dest_mm(1,:)=-XYZ_dest_mm(1,:);
                XYZ_dest_mm(2,:)=-XYZ_dest_mm(2,:);

                %  make sure coors is in 4*N size (for further transformation)
                XYZ_dest_mm = [XYZ_dest_mm; ones(1,size(XYZ_dest_mm, 2))];

        elseif contains(normMethod, 'FLIRT') || contains(normMethod, 'FSL') % Registration done by FSL (ea_flirt)
                % fuzzy match, transform can be specified as *_flirt.mat or simply *.mat
                if ~isfile(transform) % *.mat not specified
                    transform = regexprep(transform, '\.mat$', '-flirt.mat');
                    if ~isfile(transform)
                        error('Transformation file not detected!');
                    end
                end

                % vox to mm, img2imgcoord takes mm coords as input
                XYZ_src_mm = ea_get_affine(src)*XYZ_src_vx;

                % apply transform, need transpose because FSL prefer N*3 like row vector
                XYZ_dest_mm = ea_fsl_img2imgcoord(XYZ_src_mm(1:3,:)', src, dest, transform, 'linear')';

                %  make sure coors is in 4*N size (for further transformation)
                XYZ_dest_mm = [XYZ_dest_mm; ones(1,size(XYZ_dest_mm, 2))];

        else % actual FALLBACK
                % If no transformmethod is supplied, do a just-in-time
                % spm_coreg here as fallback, since linear transformation
                % is not computational expensive anyway, the transformation
                % will be saved as {$src}2{$dest}.mat
                fprintf('\nTransformation method not found!\nFallback to spm_coreg...\n');
                x = spm_coreg(dest, src);
                XYZ_dest_mm = spm_matrix(x(:)')\ea_get_affine(src)*XYZ_src_vx;

                % directory = fileparts(GetFullPath(src));
                % save([directory, filesep, ea_getmodality(src), '2', ea_getmodality(dest), '.mat'], 'x');
        end

   	% 'forwardTransform' or 'inverseTransform' stub supplied, LeadDBS's
   	% non-linear normalization, proper tranformation files will be
    % automatically detected for ANTs, FSL or SPM. 'inverseTransform'
    % automatically if you want to map src coords to dest coords
    % (internally, for ANTs, FSL and SPM, the inverse of the deformation
    % field is used for themapping).
    elseif endsWith(transform, 'forwardTransform') || endsWith(transform, 'inverseTransform')
        % Get normalization method
    	subjDir = fileparts(transform);
        options = ea_getptopts(subjDir);
        json = loadjson(options.subj.norm.log.method);
        normMethod = upper(json.method);

        if contains(normMethod, {'ANTS', 'EASYREG', 'SYNTHMORPH', 'SPM'})
            % Convert SPM deformation field to ITK format when necessary
            if contains(normMethod, 'SPM')
                ea_convert_spm_warps(options.subj);
            end

            if endsWith(transform, 'inverseTransform')
                useinverse = 1;
            else
                useinverse = 0;
            end

            % vox to mm, ANTs takes mm coords as input
            XYZ_src_mm = ea_get_affine(src)*XYZ_src_vx;

            % RAS to LPS, ANTs (ITK) use LPS coords
            XYZ_src_mm(1,:)=-XYZ_src_mm(1,:);
            XYZ_src_mm(2,:)=-XYZ_src_mm(2,:);

            % apply transform, need transpose becuase ANTs prefer N*3
            % like row vector
            try
                XYZ_dest_mm = ea_antspy_apply_transforms_to_points(XYZ_src_mm(1:3,:)', subjDir, useinverse)';
            catch
                ea_cprintf('CmdWinWarnings', 'Failed to run ANTsPy! Fallback to ANTs...\n');
                XYZ_dest_mm = ea_ants_apply_transforms_to_points(XYZ_src_mm(1:3,:)', subjDir, useinverse)';
            end

            % LPS to RAS, restore to RAS coords
            XYZ_dest_mm(1,:)=-XYZ_dest_mm(1,:);
            XYZ_dest_mm(2,:)=-XYZ_dest_mm(2,:);

            %  make sure coors is in 4*N size (for further transformation)
            XYZ_dest_mm = [XYZ_dest_mm; ones(1,size(XYZ_dest_mm, 2))];

        elseif contains(normMethod, 'FNIRT') || contains(normMethod, 'FSL')
            % if 'inverseTransform' is specified, it means
            % mapping from src coords to dest coords, i.e., NOT the
            % inverse mapping (from dest coords to src coords).
            if endsWith(transform, 'inverseTransform')
                inversemap = 0;
            else
                inversemap = 1;
            end

            % vox to mm, img2imgcoord takes mm coords as input
            XYZ_src_mm = ea_get_affine(src)*XYZ_src_vx;

            % apply transform, need transpose because FSL prefer N*3
            % like row vector
            XYZ_dest_mm = ea_fsl_apply_normalization_to_points(subjDir,XYZ_src_mm(1:3,:)',inversemap)';

            %  make sure coors is in 4*N size (for further transformation)
            XYZ_dest_mm = [XYZ_dest_mm; ones(1,size(XYZ_dest_mm, 2))];

        % % Default use SPM to do the mapping
        % elseif contains(normMethod, 'SPM')
        %     if endsWith(transform, 'inverseTransform')
        %         transform = [options.subj.norm.transform.inverseBaseName, 'spm.nii'];
        %     else
        %         transform = [options.subj.norm.transform.forwardBaseName, 'spm.nii'];
        %     end
        %
        %     XYZ_dest_mm = srcvx2destmm_deform(XYZ_src_vx, transform);
        else
            error('Normalization method not recognizable!');
        end

    % 'y_*.nii' or 'iy_*.nii' from SPM supplied, NOLINEAR case
    elseif ~isempty(regexp(transform, ['(?:\', filesep, '|^)(y|iy)_.+\.nii$'], 'once'))
        XYZ_dest_mm = srcvx2destmm_deform(XYZ_src_vx, transform);

    % '*.nii', '*.nii.gz' or '*.h5' files from ANTs, FSL or SPM (saved in ITK format) supplied, NOLINEAR case
    elseif ~isempty(regexp(transform, '\.nii$', 'once')) || ... % ANTs or FSL naming
           ~isempty(regexp(transform, '\.nii.gz$', 'once')) || ... % ANTs or FSL naming
           ~isempty(regexp(transform, '\.h5$', 'once')) % ANTs naming

        % Need to specify the transformation type
        if nargin >= 5
            normMethod = upper(varargin{5});
        else
            error('Please specify the transformation type');
        end

        switch normMethod

            case {'ANTS', 'EASYREG', 'SYNTHMORPH', 'SPM'} % ANTs or SPM used
                if nargin >= 6
                    useinverse = varargin{6};
                else
                    useinverse = 0; % suppose proper deformation field specified, no need to invert
                end

                % vox to mm, ANTs takes mm coords as input
                XYZ_src_mm = ea_get_affine(src)*XYZ_src_vx;

                % RAS to LPS, ANTs (ITK) use LPS coords
                XYZ_src_mm(1,:)=-XYZ_src_mm(1,:);
                XYZ_src_mm(2,:)=-XYZ_src_mm(2,:);

                % apply transform, need transpose becuase ANTs prefer N*3
                % like row vector
                try
                    XYZ_dest_mm = ea_antspy_apply_transforms_to_points(XYZ_src_mm(1:3,:)', transform, useinverse)';
                catch
                    ea_cprintf('CmdWinWarnings', 'Failed to run ANTsPy! Fallback to ANTs...\n');
                    XYZ_dest_mm = ea_ants_apply_transforms_to_points(XYZ_src_mm(1:3,:)', transform, useinverse)';
                end

                % LPS to RAS, restore to RAS coords
                XYZ_dest_mm(1,:)=-XYZ_dest_mm(1,:);
                XYZ_dest_mm(2,:)=-XYZ_dest_mm(2,:);

                %  make sure coors is in 4*N size (for further transformation)
                XYZ_dest_mm = [XYZ_dest_mm; ones(1,size(XYZ_dest_mm, 2))];

            case {'FSL', 'FNIRT'} % FSL (ea_fnirt) used
                % vox to mm, img2imgcoord takes mm coords as input
                XYZ_src_mm = ea_get_affine(src)*XYZ_src_vx;

                % apply transform, need transpose because FSL prefer N*3
                % like row vector
                XYZ_dest_mm = ea_fsl_img2imgcoord(XYZ_src_mm(1:3,:)', src, dest, transform, 'nonlinear')';

                %  make sure coors is in 4*N size (for further transformation)
                XYZ_dest_mm = [XYZ_dest_mm; ones(1,size(XYZ_dest_mm, 2))];

            otherwise
                error(['Unsupported transformation type: ', normMethod])
        end

    else
        error(['Unsupported transformation file:\n', transform])
    end

    % Optional output from dest vx coords
    if nargout == 2
        if nargin >= 4 % dest image specified
            XYZ_dest_vx = ea_get_affine(dest) \ XYZ_dest_mm;

        elseif ~isempty(regexp(transform, 'sn\.mat$', 'once')) % transform is SPM DCT file
            sn = load(transform, 'VF');
            XYZ_dest_vx = sn.VF.mat \ XYZ_dest_mm;

        elseif ~isempty(regexp(transform, '_spm\.mat$', 'once')) % transform is from ea_spm_coreg
            affine = load(transform, 'fixed');
            XYZ_dest_vx = affine.fixed \ XYZ_dest_mm;
        end

        XYZ_dest_vx=XYZ_dest_vx(1:3,:);
    end

    XYZ_dest_mm=XYZ_dest_mm(1:3,:);
end


function coord = srcvx2destmm_sn(coord, matname)
% returns mm coordinates based on the old version of the SPM tranformation:
% '*_sn.mat'

sn = load(matname);
Tr = sn.Tr;

if numel(Tr) ~= 0 % DCT warp: src_vox displacement
    d = sn.VG(1).dim(1:3); % (since VG may be 3-vector of TPM volumes)
    dTr = size(Tr);
    basX = spm_dctmtx(d(1), dTr(1), coord(1,:)-1);
    basY = spm_dctmtx(d(2), dTr(2), coord(2,:)-1);
    basZ = spm_dctmtx(d(3), dTr(3), coord(3,:)-1);
    for i = 1:size(coord, 2)
        bx = basX(i, :);
        by = basY(i, :);
        bz = basZ(i, :);
        tx = reshape(...
            reshape(Tr(:,:,:,1),dTr(1)*dTr(2),dTr(3))*bz',dTr(1),dTr(2) );
        ty = reshape(...
            reshape(Tr(:,:,:,2),dTr(1)*dTr(2),dTr(3))*bz',dTr(1),dTr(2) );
        tz =  reshape(...
            reshape(Tr(:,:,:,3),dTr(1)*dTr(2),dTr(3))*bz',dTr(1),dTr(2) );
        coord(1:3,i) = coord(1:3,i) + [bx*tx*by' ; bx*ty*by' ; bx*tz*by'];
    end
end

% Affine: src_vx (possibly displaced by above DCT) to dest_vx
coord = sn.VF.mat * sn.Affine * coord;


function dest_mm = srcvx2destmm_deform(src_vx, deform)
% returns mm coordinates based on deformation field file 'y_*.nii' from src
% image to dest image

if ischar(deform)
    deform = spm_vol([repmat(deform,3,1),[',1,1';',1,2';',1,3']]);
end

src_vx = double(src_vx);
dest_mm = [spm_sample_vol(deform(1,:),src_vx(1,:),src_vx(2,:),src_vx(3,:),1);...
          spm_sample_vol(deform(2,:),src_vx(1,:),src_vx(2,:),src_vx(3,:),1);...
          spm_sample_vol(deform(3,:),src_vx(1,:),src_vx(2,:),src_vx(3,:),1)];
if size(src_vx,1) == 4
    dest_mm = [dest_mm; src_vx(4,:)];
end
