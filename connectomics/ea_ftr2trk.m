function ea_ftr2trk(ftrfile, specs, LPS)
% export FTR matrix to TrackVis trk format
%
% specs can also be the path of the nifti file which defines the space. If
% not specified or empty, MNI space T1 will be used as reference.
%
% If the trk is going to be visualized in Surf-Ice, LPS should be set to 1
% to fix the orientation.

[directory, ftrname, ext] = fileparts(ftrfile);
if isempty(directory)
    directory = '.';
end

if isempty(ext)
    ftrfile = [ftrfile, '.mat'];
end

disp('Loading FTR-File...');
[fibs, idx, voxmm] = ea_loadfibertracts(ftrfile);

% Convert ONE-BASED indexing to ZERO-BASED indexing
if strcmp(voxmm,'vox')
    fibs(:,1:3) = fibs(:,1:3) - 1;
end

%% set header
header.id_string = ['TRACK', char(0)];
header.voxel_order = ['LPS', char(0)];
header.pad1 = repmat(char(0), 1, 2);

if ~exist('specs','var') || isempty(specs) % Use MNI T1 as reference space by default.
    disp('Header from space template ...');
    spacedef = ea_getspacedef;
    refhdr = ea_fslhd([ea_space, spacedef.templates{1}, '.nii']);
    specs.origin = [0,0,0];
    specs.dim = [refhdr.dim1, refhdr.dim2, refhdr.dim3];
    specs.affine = ea_get_affine([ea_space, spacedef.templates{1}, '.nii'], 0);
    header.pad2 = ['RAS', char(0)];
elseif isstruct(specs)
    % Suppose that the affine matrix is from SPM
    specs.affine(:,4) = specs.affine(:,4) + sum(specs.affine(:,1:3),2);
    header.pad2 = [ea_aff2axcodes(specs.affine), char(0)];
elseif isfile(specs) % Use the specified nifti as reference space.
    disp(['Header from ',specs,' ...']);
    refimage = specs;
    refhdr = ea_fslhd(refimage);
    specs = struct;
    specs.origin = [0,0,0];
    specs.dim = [refhdr.dim1, refhdr.dim2, refhdr.dim3];
    specs.affine = ea_get_affine(refimage, 0);
    header.pad2 = [ea_aff2axcodes(specs.affine), char(0)];
end

specs = ea_aff2hdr(specs.affine, specs, 1);

header.dim = specs.dim;
header.voxel_size = specs.voxel_size;
header.vox_to_ras  =  specs.vox_to_ras;
header.image_orientation_patient = specs.image_orientation_patient;
header.origin = [0 0 0]; % as doc says, trackvis will always use 0 0 0 as origin.

header.n_scalars = 0;
header.scalar_name = char(repmat(' ',10,20));
header.n_properties = 0;
header.property_name = char(repmat(' ',10,20));
header.reserved = char(repmat(' ',444,1));
header.invert_x = 0;
header.invert_y = 0;
header.invert_z = 0;
header.swap_xy = 0;
header.swap_yz = 0;
header.swap_zx = 0;
header.n_count = length(idx);% header.invert_x = 1;
header.version = 2;
header.hdr_size = 1000;

%% convert data
disp('Constructing data...');
tracks = struct('nPoints',nan,'matrix',nan);
offset = 1;
for track_number=1:length(idx)
    tracks(1,track_number).nPoints = idx(track_number);
    tracks(1,track_number).matrix = fibs(offset:offset+idx(track_number)-1,1:3);
    offset = offset+idx(track_number);
end

if strcmp(voxmm,'mm') % have to retranspose to vox
    disp('mm to vox conversion...');
    for i=1:length(tracks)
        tracks(i).matrix = [tracks(i).matrix,ones(size(tracks(i).matrix,1),1)]';
        tracks(i).matrix = specs.affine\tracks(i).matrix;
        if exist('LPS', 'var') && LPS
            tracks(i).matrix(1,:) = refhdr.dim1-1-tracks(i).matrix(1,:);
            tracks(i).matrix(2,:) = refhdr.dim2-1-tracks(i).matrix(2,:);
        end
        tracks(i).matrix = tracks(i).matrix(1:3,:)';
    end
end

for i = 1:length(tracks)
    try
        tracks(i).matrix = bsxfun(@times, tracks(i).matrix,header.voxel_size);
    catch
        tracks(i).matrix = bsxfun(@times, tracks(i).matrix',header.voxel_size);
    end
end

%% write .trk file
disp('Writing trk file...');
ea_trk_write(header,tracks,[directory,filesep,ftrname,'.trk']);

disp('Conversion finished.');


function trk_hdr = ea_aff2hdr(affine, trk_hdr, pos_vox, set_order)
% Set affine matrix into trackvis header 'trk_hdr'
%
% pos_vos : None or bool
%     If None, currently defaults to False. If False, allow negative voxel
%     sizes in header to record axis flips. Negative voxels cause problems
%     for TrackVis.  If True, enforce positive voxel sizes.
% set_order : None or bool
%     If None, currently defaults to False. If False, do not set 'voxel_order'
%     field in 'trk_hdr'. If True, calculcate 'voxel_order' from 'affine' and
%     set into 'trk_hdr'.
%
% Notes
% -----
% version 2 of the trackvis header has a dedicated field for the nifti RAS
% affine. In theory trackvis 1 has enough information to store an affine, with
% the fields 'origin', 'voxel_size' and 'image_orientation_patient'.
% Unfortunately, to be able to store any affine, we'd need to be able to set
% negative voxel sizes, to encode axis flips. This is because
% 'image_orientation_patient' is only two columns of the 3x3 rotation matrix,
% and we need to know the number of flips to reconstruct the third column
% reliably. It turns out that negative flips upset trackvis (the
% application). The application also ignores the origin field, and may not
% use the 'image_orientation_patient' field.
%

if nargin < 3
    pos_vox = false;
end

if nargin < 4
    set_order = false;
end

try
    version = trk_hdr.version;
catch
    version = 2;
end

if version == 2
    vox_to_ras = affine;
    % Check orientation
    % Adapt the affine matrix if it's not in RAS orientation.
    if vox_to_ras(1) < 0
        vox_to_ras(1,:) = vox_to_ras(1,:) * -1;
    end
    if vox_to_ras(6) < 0
        vox_to_ras(2,:) = vox_to_ras(2,:) * -1;
    end
    if vox_to_ras(11) < 0
        vox_to_ras(3,:) = vox_to_ras(3,:) * -1;
    end
    trk_hdr.vox_to_ras = vox_to_ras;
end

if set_order
    trk_hdr.voxel_order = ea_aff2axcodes(affine);
end

% affine to go from DICOM LPS to MNI RAS space
DPCS_TO_TAL = diag([-1,-1,1,1]);
% Now on dodgy ground with DICOM fields in header
% RAS to DPCS output
affine = DPCS_TO_TAL*affine;
trans = affine(1:3,4);
% Get zooms
RZS = affine(1:3, 1:3);
zooms = sqrt(sum(RZS.*RZS));
RS = RZS./repmat(zooms,3,1);
% If you said we could, adjust zooms to make RS correspond (below) to a true
% rotation matrix.  We need to set the sign of one of the zooms to deal with
% this. TrackVis doesn't like negative zooms at all, so you might want to
% disallow this with the pos_vox option.
if ~pos_vox && det(RS)<0
    zooms(1) = zooms(1)*-1;
    RS(:,1) = RS(:,1)*-1;
end
% retrieve rotation matrix from RS with polar decomposition.
% Discard shears because we cannot store them.
[P, ~, Qs] = svd(RS);
R = P*Qs';
% check if R is an orthogonal matrix
assert(ea_allclose(R*R',eye(3)), 'non-orthogonal R matrix')

trk_hdr.origin = trans;
trk_hdr.voxel_size = zooms;
trk_hdr.image_orientation_patient = reshape(R(:,1:2),1,[]);
