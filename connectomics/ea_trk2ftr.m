function [fibers,idx] = ea_trk2ftr(trk_in,type)
% type can be:
% 1 - 'DSI Studio / QSDR'
% 2 - 'Normative Connectome / 2009b space'
% 3 - 'Select .nii file'
% or a filepath (string) pointing to the .nii file associated with the
% tracts.
% in case of 1 (DSI Studio / QSDR), heuristic transforms specified under
% http://dsi-studio.labsolver.org/Manual/Reconstruction#TOC-Q-Space-Diffeomorphic-Reconstruction-QSDR-
% are applied. In all other cases, vox2mm transform of a .nii header will
% be applied.

% transforms .trk to fibers.
% CAVE: fiber information in the trk needs to conform to the MNI space!
if ~exist('type', 'var')
    type = nan;
end

%% read .trk file
[~, fn, ext] = fileparts(trk_in);
if strcmp(ext,'.gz')
    uuid = ea_generate_uuid;
    td = [ea_getleadtempdir, uuid];
    mkdir(td);
    gunzip(trk_in, td);
    trk_in = fullfile(td, fn);
    [header,tracks] = ea_trk_read(trk_in);
    rmdir(td, 's');
else
    [header, tracks] = ea_trk_read(trk_in);
end

%% create variables needed
ea_fibformat = '1.0';
idx = [tracks(:).nPoints]';
idx2 = cumsum(idx);
fibers = NaN(4, sum(idx));

start = 1;
fibercount = numel(idx);
ea_dispercent(0, 'Converting trk to fibers.');
for a=1:fibercount
    ea_dispercent(a/fibercount);
    stop = idx2(a);
    fibers(1:3, start:stop) = tracks(a).matrix';
    fibers(4, start:stop) = a;
    start = stop+1;
end
ea_dispercent(1,'end');

fibers(1:3,:) = fibers(1:3,:)./header.voxel_size';

%% transform fibers to template origin
if isnan(type)
    answ = questdlg('Please select type of .trk data', 'Choose origin of .trk data', 'DSI Studio / QSDR', 'Normative Connectome / 2009b space', 'Select .nii file', 'DSI Studio / QSDR');
else
    if ischar(type)
        answ = 'specified_image';
    else
        if type==1
            answ = 'DSI Studio / QSDR';
        elseif type==2
            answ = 'Normative Connectome / 2009b space';
        elseif type==3
            answ = 'Select .nii file';
        end
    end
end

switch answ
    case 'DSI Studio / QSDR'
        % http://dsi-studio.labsolver.org/Manual/Reconstruction#TOC-Q-Space-Diffeomorphic-Reconstruction-QSDR-
        fibers(1,:) = 78.0 - fibers(1,:);
        fibers(2,:) = 76.0 - fibers(2,:);
        fibers(3,:) = -50.0 + fibers(3,:);
    otherwise
        switch answ
            case 'specified_image'
                nii = ea_load_nii(type); % load in image supplied to function.
            case 'Normative Connectome / 2009b space'
                spacedef = ea_getspacedef;
                nii = ea_load_nii([ea_space, spacedef.templates{1}]);
            case 'Select .nii file'
                [fname, pathname] = uigetfile({'.nii','.nii.gz'}, 'Choose Nifti file for space definition');
                nii = ea_load_nii(fullfile(pathname, fname));
        end

        tfib = [fibers(1:3,:); ones(1, size(fibers,2))];

        tmat = header.vox_to_ras;
        if isempty(find(tmat, 1))
            % method Andreas (DSIStudio .trk)
            tmat(1:3,4) = header.dim';
        else
            % method Till, in case vox_to_ras is null matrix (TrackVis?)
            tmat = nii.mat;
            tmat(1,1) = 1;
            tmat(2,2) = 1;
            tmat(3,3) = 1;
        end
        tfib = tmat*tfib;
        fibers(1:3,:) = tfib(1:3,:);
end

fibers = single(fibers');
