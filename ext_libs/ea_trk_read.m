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
if ~isfile(filePath)
    ea_error('TRK file not found!', showdlg=0, simpleStack=1);
end

fid    = fopen(filePath, 'r');
header = get_header(fid);

% Check for byte order
if header.hdr_size~=1000
    fclose(fid);
    fid    = fopen(filePath, 'r', 'b'); % Big endian for old PPCs
    header = get_header(fid);
end

if header.hdr_size~=1000
    ea_error('TRK Header length is wrong!', showdlg=0, simpleStack=1);
end

if ~any(header.image_orientation_patient)
    % ea_cprintf('CmdWinWarnings', '"image_orientation_patient" was not set properly in the header!\nWill try to set it in a heuristic way.\n');

    % if header.vox_to_ras(1) > 0
    %     header.image_orientation_patient(1:3) = [-1 0 0];
    % else
    %     header.image_orientation_patient(1:3) = [1 0 0];
    % end
    % 
    % if header.vox_to_ras(6) > 0
    %     header.image_orientation_patient(4:6) = [0 -1 0];
    % else
    %     header.image_orientation_patient(4:6) = [0 1 0];
    % end

    ea_cprintf('CmdWinWarnings', '"image_orientation_patient" was not set properly in the header!\nFallback to [1 0 0 0 1 0] (LPS).\n');
    header.image_orientation_patient = [1 0 0 0 1 0];
end

if ~any(ismember(header.pad2(1:3), 'RASLPI'))
    ea_cprintf('CmdWinWarnings', '"pad2" was not set properly in the header!\nWill try to set it in a heuristic way.\n');
    
    if header.vox_to_ras(1) > 0
        header.pad2(1) = 'R';
    else
        header.pad2(1) = 'L';
    end

    if header.vox_to_ras(6) > 0
        header.pad2(2) = 'A';
    else
        header.pad2(2) = 'P';
    end

    if header.vox_to_ras(11) > 0
        header.pad2(3) = 'S';
    else
        header.pad2(3) = 'I';
    end
end

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
        coords(:,ix) = (header.dim(ix)-1) * header.voxel_size(ix) - coords(:,ix);
    end
    if header.image_orientation_patient(3+iy) < 0
        coords(:,iy) = (header.dim(iy)-1) * header.voxel_size(iy) - coords(:,iy);
    end
    tracks(iTrk).matrix(:,1:3) = coords;
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
