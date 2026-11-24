
function [NiiFileOut, FlipB] = LeG_rotateImg2Standard(NiiFile, varargin)
% Reorient a NIfTI image to the orientation expected by LeGUI:
% voxel dimension 1 -> negative world X
% voxel dimension 2 -> positive world Y
% voxel dimension 3 -> positive world Z
%
% IMPORTANT:
% Permutes/flips voxel data and updates the affine consistently.
% Preserves the physical world position of every voxel.
% Does not modify the input image.
%
% Output: sd-prefixed NIfTI, plus indicator of X-axis reversal.

    V = spm_vol(NiiFile);
    img = spm_read_vols(V);

    oldMat = V.mat;
    A = oldMat(1:3,1:3);
    oldSize = size(img);
    oldSize(end+1:3) = 1;

    % Desired voxel directions in world coordinates
    desiredSigns = [-1, 1, 1];

    % Determine which voxel axis corresponds to each world axis.
    axisLength = sqrt(sum(A.^2,1));

    if any(~isfinite(axisLength)) || any(axisLength <= 0)
        error('Invalid input NIfTI affine.');
    end

    dirs = A ./ axisLength;

    [strength, oldAxis] = max(abs(dirs), [], 2);
    oldAxis = oldAxis(:)';

    % This method is intended for approximately axis-aligned
    % images. Oblique images require a different strategy.
    if numel(unique(oldAxis)) ~= 3 || any(strength < 0.98)
        error(['Image is oblique or has ambiguous axes. ' ...
               'Cannot safely reorient using 90-degree rotations.']);
    end

    % Reorder data so dimensions correspond to world X,Y,Z.
    imgR = permute(img, oldAxis);

    % T maps NEW voxel coordinates back to OLD voxel coordinates.
    % SPM uses 1-based voxel indexing.
    T = zeros(4);
    T(4,4) = 1;

    FlipB = false;

    for k = 1:3
        j = oldAxis(k);

        currentSign = sign(A(k,j));
        needsFlip = currentSign ~= desiredSigns(k);

        if needsFlip
            imgR = flip(imgR, k);

            T(j,k) = -1;
            T(j,4) = oldSize(j) + 1;
        else
            T(j,k) = 1;
        end

        if k == 1
            FlipB = needsFlip;
        end
    end

    % Update affine using exactly the same voxel transform.
    newMat = oldMat * T;

    % Sanity check: physical voxel positions must be preserved.
    newSize = size(imgR);
    newSize(end+1:3) = 1;

    testVox = [
        1, 1, 1;
        newSize(1), 1, 1;
        1, newSize(2), 1;
        1, 1, newSize(3);
        newSize(1:3)
    ];

    newVox = [testVox, ones(size(testVox,1),1)]';

    oldWorld = oldMat * T * newVox;
    newWorld = newMat * newVox;

    err = max(abs(oldWorld(:) - newWorld(:)));

    if err > 1e-6
        error('Reorientation validation failed.');
    end

    % Save as sd* without modifying the original.
    [folder, name, ext] = fileparts(V.fname);

    if startsWith(name, 'sd')
        error('Input already has sd prefix: %s', V.fname);
    end

    NiiFileOut = fullfile(folder, ['sd', name, ext]);

    if isfile(NiiFileOut)
        error('Output already exists: %s', NiiFileOut);
    end

    V.fname = NiiFileOut;
    V.mat = newMat;
    V.dim = newSize(1:3);

    spm_write_vol(V, imgR);

    fprintf('\nReorientation completed.\n');
    fprintf('X-axis flip applied: %d\n', FlipB);
    fprintf('Voxel/world coordinate error: %.3g mm\n', err);
    fprintf('Output: %s\n', NiiFileOut);

end
