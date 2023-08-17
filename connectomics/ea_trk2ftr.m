function [fibers, idx] = ea_trk2ftr(trkFile, ref, outputFile)
% Convert trk to ftr (fibers format in Lead-DBS)
%
% ref can be:
%   1        - Lead-DBS / MNI152 NLin 2009b Asym
%   2        - Choose a reference NIfTI file
%   3        - Interactively select ref type (1 or 2 above)
%   [string] - Path of a NIfTI file defining the reference space

if ~exist('ref', 'var')
    ref = 'mni';
end

if isnumeric(ref)
    ref = num2str(ref);
end

if ~exist('outputFile', 'var')
    outputFile = 0;
end

% read .trk file
[~, fn, ext] = fileparts(trkFile);
if strcmp(ext,'.gz')
    uuid = ea_generate_uuid;
    td = [ea_getleadtempdir, uuid];
    mkdir(td);
    gunzip(trkFile, td);
    trkFile = fullfile(td, fn);
    [header,tracks] = ea_trk_read(trkFile);
    rmdir(td, 's');
else
    [header, tracks] = ea_trk_read(trkFile);
end

% Get affine from reference
if ismember(ref, {'3', 'select', 'ask', 'interactive'})
    answer = questdlg('Please specify the reference space of the trk file.', '', 'Lead-DBS / MNI152 NLin 2009b Asym', 'Choose a reference NIfTI file', 'Lead-DBS / MNI152 NLin 2009b Asym');
    if isempty(answer)
        return;
    elseif strcmp(answer, 'Lead-DBS / MNI152 NLin 2009b Asym')
        ref = 'mni';
    elseif strcmp(answer, 'Choose a reference NIfTI file')
        ref = 'nii';
    end
end

switch ref
    case {'1', 'lead', 'leaddbs', 'mni'}
        % MNI space trk generated within Lead-DBS already has the vox_to_ras properly set
        affine = header.vox_to_ras;
    case {'2', 'nii', 'nifti', 'other'}
        [fname, pathname] = uigetfile({'*.nii','*.nii.gz'}, 'Choose a reference NIfTI file');
        if ischar(fname)
            affine = ea_get_affine(fullfile(pathname, fname), 0);
        else
            return;
        end
    otherwise
        if isfile(ref)
            affine = ea_get_affine(ref, 0);
        else
            ea_error('Unrecognized reference type!', 'simpleStack', true);
        end
end

% Contruct fibers and idx needed for ftr mat
idx = vertcat(tracks.nPoints);
fibers = vertcat(tracks.matrix);
fibers(:, 4) = repelem(1:length(idx), idx)';

% Rescale
fibers(:, 1:3) = fibers(:, 1:3)./header.voxel_size;

% Flip L/R, P/A and I/S in case affine and pad2 doesn't match 
if strcmp(header.pad2(1), 'L') && affine(1) > 0 || strcmp(header.pad2(1), 'R') && affine(1) < 0
    % Flip L and R
    fibers(:, 1) = header.dim(1) - 1 - fibers(:, 1);
end

if strcmp(header.pad2(2), 'P') && affine(6) > 0 || strcmp(header.pad2(2), 'A') && affine(6) < 0
    % Flip P and A
    fibers(:, 2) = header.dim(2) - 1 - fibers(:, 2);
end

if strcmp(header.pad2(3), 'I') && affine(11) > 0 || strcmp(header.pad2(3), 'S') && affine(11) < 0
    % Flip I and S
    fibers(:, 3) = header.dim(3) - 1 - fibers(:, 3);
end

% Convert voxel to mm 
fibers(:, 1:3) = ea_vox2mm(fibers(:, 1:3), affine);

fibers = single(fibers);

% Optionally save ftr mat
if outputFile
    outputFile = replace(erase(trkFile, '.gz'), '.trk', '.mat');
    if isfile(outputFile)
        answer = questdlg('File already exists!', '', 'Overwrite', 'Specify a New Name', 'Overwrite');
        if isempty(answer)
            return;
        elseif strcmp(answer, 'Specify a New Name')
            [fname, pathname] = uiputfile({'*.mat'}, 'Choose a reference NIfTI file', outputFile);
            outputFile = fullfile(pathname, fname);
        end
    end

    ea_fibformat = '1.0';
    fourindex = 1;
    voxmm = 'mm';
    save(outputFile, 'ea_fibformat', 'fibers', 'fourindex', 'idx', 'voxmm', '-v7.3');
end
