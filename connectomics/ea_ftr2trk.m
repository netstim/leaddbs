function ea_ftr2trk(ftrfilename,directory,specs)
% export FTR matrix to TrackVis trk format

if ischar(ftrfilename)
    disp('Loading FTR-File...');
    [fibs,idx,voxmm,mat]=ea_loadfibertracts([directory,ftrfilename,'.mat']);
else % direct ftr import
    ea_error('Direct FTR import not supported at present.');
    fibs=ftrfilename{1};
    ftrfilename{1}=[];
end

%% set header
[header, ~]=ea_trk_read([ea_getearoot,'ext_libs',filesep,'example.trk']);
if strcmp(voxmm,'vox')
    if ~isempty(mat)
        if isempty(specs.affine)
            specs.affine=mat;
        else
            if ~isequal(mat,specs.affine)
                ea_error('Affine matrix of Fibertracts and image do not match');
            end
        end
    end
end

% check if x-axis of the affine matrix is negative, flip it if so
if det(specs.affine) < 0
    specs.affine = diag([-1 1 1 1])*specs.affine;
end
specs = ea_aff2hdr(specs.affine, specs);

try
    header.voxel_size=specs.voxel_size;
catch
    header.voxel_size=fs.vox;
end
header.dim=specs.dim;
header.origin=[0 0 0]; % as doc says, trackvis will always use 0 0 0 as origin.
header.n_scalars=0;
header.scalar_name=char(repmat(' ',10,20));
header.n_properties=0;
header.property_name=char(repmat(' ',10,20));
header.reserved=char(repmat(' ',444,1));
header.image_orientation_patient=specs.image_orientation_patient;
header.invert_x=0;
header.invert_y=0;
header.invert_z=0;
header.swap_xy=0;
header.swap_yz=0;
header.swap_zx=0;
header.n_count=length(idx);% header.invert_x=1;
header.version=2;
header.hdr_size=1000;

%% convert data
tracks=struct('nPoints',nan,'matrix',nan);
offset=1;
for track_number=1:length(idx)
   tracks(1,track_number).nPoints=idx(track_number);
   tracks(1,track_number).matrix=fibs(offset:offset+idx(track_number)-1,1:3);
   offset=offset+idx(track_number);
end

if strcmp(voxmm,'mm') % have to retranspose to vox
    for i=1:length(tracks)
       tracks(i).matrix=[tracks(i).matrix,ones(size(tracks(i).matrix,1),1)]';
        tracks(i).matrix=specs.affine\tracks(i).matrix;
        tracks(i).matrix=tracks(i).matrix(1:3,:)';
    end
end

for i = 1:length(tracks)
    try
        tracks(i).matrix=bsxfun(@times, tracks(i).matrix,header.voxel_size);
    catch
        tracks(i).matrix=bsxfun(@times, tracks(i).matrix',header.voxel_size);
    end
end

%% write .trk file
if ischar(ftrfilename)
    ea_trk_write(header,tracks,[directory,ftrfilename,'.trk']);
else
    ea_trk_write(header,tracks,[directory,ftrfilename{2},'.trk']);
end


function [header,tracks] = ea_trk_read(filePath)
%TRK_READ - Load TrackVis .trk files
%TrackVis displays and saves .trk files in LPS orientation. After import, this
%function attempts to reorient the fibers to match the orientation of the
%original volume data.
%
% Syntax: [header,tracks] = trk_read(filePath)
%
% Inputs:
%    filePath - Full path to .trk file [char]
%
% Outputs:
%    header - Header information from .trk file [struc]
%    tracks - Track data structure array [1 x nTracks]
%      nPoints - # of points in each streamline
%      matrix  - XYZ coordinates and associated scalars [nPoints x 3+nScalars]
%      props   - Properties of the whole tract (ex: length)
%
% Example:
%    exDir           = '/path/to/along-tract-stats/example';
%    subDir          = fullfile(exDir, 'subject1');
%    trkPath         = fullfile(subDir, 'CST_L.trk');
%    [header tracks] = trk_read(trkPath);
%
% Other m-files required: none
% Subfunctions: get_header
% MAT-files required: none
%
% See also: http://www.trackvis.org/docs/?subsect=fileformat
%           http://github.com/johncolby/along-tract-stats/wiki/orientation

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Mar 2010

% Parse in header
fid    = fopen(filePath, 'r');
header = get_header(fid);

% Check for byte order
if header.hdr_size~=1000
    fclose(fid);
    fid    = fopen(filePath, 'r', 'b'); % Big endian for old PPCs
    header = get_header(fid);
end

if header.hdr_size~=1000, ea_error('FTR-Header length is wrong'), end

% Check orientation
[~, ix] = max(abs(header.image_orientation_patient(1:3)));
[~, iy] = max(abs(header.image_orientation_patient(4:6)));
iz = 1:3;
iz([ix iy]) = [];

% Parse in body
tracks(header.n_count).nPoints = 0;

for iTrk = 1:header.n_count
    tracks(iTrk).nPoints = fread(fid, 1, 'int');
    tracks(iTrk).matrix  = fread(fid, [3+header.n_scalars, tracks(iTrk).nPoints], 'float')';
    if header.n_properties
        tracks(iTrk).props = fread(fid, header.n_properties, 'float');
    end

    % Modify orientation of tracks (always LPS) to match orientation of volume
    header.dim        = header.dim([ix iy iz]);
    header.voxel_size = header.voxel_size([ix iy iz]);
    coords = tracks(iTrk).matrix(:,1:3);
    coords = coords(:,[ix iy iz]);
    if header.image_orientation_patient(ix) < 0
        coords(:,ix) = header.dim(ix)*header.voxel_size(ix) - coords(:,ix);
    end
    if header.image_orientation_patient(3+iy) < 0
        coords(:,iy) = header.dim(iy)*header.voxel_size(iy) - coords(:,iy);
    end
    tracks(iTrk).matrix(:,1:3) = coords;
end

fclose(fid);


function header = get_header(fid)

header.id_string                 = fread(fid, 6, '*char')';
header.dim                       = fread(fid, 3, 'short')';
header.voxel_size                = fread(fid, 3, 'float')';
header.origin                    = fread(fid, 3, 'float')';
header.n_scalars                 = fread(fid, 1, 'short')';
header.scalar_name               = fread(fid, [20,10], '*char')';
header.n_properties              = fread(fid, 1, 'short')';
header.property_name             = fread(fid, [20,10], '*char')';
header.vox_to_ras                = fread(fid, [4,4], 'float')';
header.reserved                  = fread(fid, 444, '*char');
header.voxel_order               = fread(fid, 4, '*char')';
header.pad2                      = fread(fid, 4, '*char')';
header.image_orientation_patient = fread(fid, 6, 'float')';
header.pad1                      = fread(fid, 2, '*char')';
header.invert_x                  = fread(fid, 1, 'uchar');
header.invert_y                  = fread(fid, 1, 'uchar');
header.invert_z                  = fread(fid, 1, 'uchar');
header.swap_xy                   = fread(fid, 1, 'uchar');
header.swap_yz                   = fread(fid, 1, 'uchar');
header.swap_zx                   = fread(fid, 1, 'uchar');
header.n_count                   = fread(fid, 1, 'int')';
header.version                   = fread(fid, 1, 'int')';
header.hdr_size                  = fread(fid, 1, 'int')';


function ea_trk_write(header,tracks,savePath)
%TRK_WRITE - Write TrackVis .trk files
%
% Syntax: trk_write(header,tracks,savePath)
%
% Inputs:
%    header   - Header information for .trk file [struc]
%    tracks   - Track data struc array [1 x nTracks]
%      nPoints  - # of points in each track
%      matrix   - XYZ coordinates and associated scalars [nPoints x 3+nScalars]
%      props    - Properties of the whole tract
%    savePath - Path where .trk file will be saved [char]
%
% Output files:
%    Saves .trk file to disk at location given by 'savePath'.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRK_READ

% Author: John Colby (johncolby@ucla.edu)
% UCLA Developmental Cognitive Neuroimaging Group (Sowell Lab)
% Apr 2010

fid = fopen(savePath, 'w');

% Write header
fwrite(fid, header.id_string, '*char');
fwrite(fid, header.dim, 'short');
fwrite(fid, header.voxel_size, 'float');
fwrite(fid, header.origin, 'float');
fwrite(fid, header.n_scalars , 'short');
fwrite(fid, header.scalar_name', '*char');
fwrite(fid, header.n_properties, 'short');
fwrite(fid, header.property_name', '*char');
fwrite(fid, header.vox_to_ras', 'float');
fwrite(fid, header.reserved, '*char');
fwrite(fid, header.voxel_order, '*char');
fwrite(fid, header.pad2, '*char');
fwrite(fid, header.image_orientation_patient, 'float');
fwrite(fid, header.pad1, '*char');
fwrite(fid, header.invert_x, 'uchar');
fwrite(fid, header.invert_y, 'uchar');
fwrite(fid, header.invert_z, 'uchar');
fwrite(fid, header.swap_xy, 'uchar');
fwrite(fid, header.swap_yz, 'uchar');
fwrite(fid, header.swap_zx, 'uchar');
fwrite(fid, header.n_count, 'int');
fwrite(fid, header.version, 'int');
fwrite(fid, header.hdr_size, 'int');

% Check orientation
[tmp ix] = max(abs(header.image_orientation_patient(1:3)));
[tmp iy] = max(abs(header.image_orientation_patient(4:6)));
iz = 1:3;
iz([ix iy]) = [];

% Write body
for iTrk = 1:header.n_count
    % Modify orientation back to LPS for display in TrackVis
    header.dim        = header.dim([ix iy iz]);
    header.voxel_size = header.voxel_size([ix iy iz]);
    coords = tracks(iTrk).matrix(:,1:3);
    coords = coords(:,[ix iy iz]);
    if header.image_orientation_patient(ix) < 0
        coords(:,ix) = header.dim(ix)*header.voxel_size(ix) - coords(:,ix);
    end
    if header.image_orientation_patient(3+iy) < 0
        coords(:,iy) = header.dim(iy)*header.voxel_size(iy) - coords(:,iy);
    end
    tracks(iTrk).matrix(:,1:3) = coords;

    fwrite(fid, tracks(iTrk).nPoints, 'int');
    fwrite(fid, tracks(iTrk).matrix', 'float');
    if header.n_properties
        fwrite(fid, tracks(iTrk).props, 'float');
    end
end

fclose(fid);


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
    trk_hdr.vox_to_ras = affine;
end

if set_order
    trk_hdr.voxel_order = aff2axcodes(affine);
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


function close = ea_allclose(a, b, rtol, atol)
% Determine if two arrays are element-wise equal within a tolerance.

if nargin < 3
    rtol = 1e-05;
end
if nargin < 4
    atol = 1e-08;
end

close = all( abs(a(:)-b(:)) <= atol+rtol*abs(b(:)) );

