% Simplified version of flirtmat2worldmat removing the various dependencies
% to SPM. Instead using PaCER's Nifti Classes which in turn rely on the
% Nifti Toolbox. Thus it works with .nii.gz files too and has a lot less
% dependencies in general.
% 2017 Andreas Husch.
%
% Based on  flirtmat2worldmat by Ged Rigway. See Original
% Text below:
%
% flirtmat2worldmat: convert saved flirt matrix to NIfTI world coords matrix
% flirt matrix is from text file specified in "flirt -omat mat.txt" command
% world matrix maps from NIfTI world coordinates in target to source. Note:
% mat.txt contains a mapping from source to target in FSL's *scaled* coords
% which are not NIfTI world coordinates, and note src-trg directionality!
% worldmat from this script reproduces "img2imgcoord -mm ...".
%
% The script can also return a matrix to map from target to source voxels
% in MATLAB/SPM's one-based convention, or in FSL's zero-based convention
%
% Example:
%  [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmat(flirtmat, src, trg);
%
% See also: worldmat2flirtmat, flirtmat_read, flirtmat_write

% Copyright 2009 Ged Ridgway <ged.ridgway gmail.com>

function [worldmat spmvoxmat fslvoxmat] = flirtmat2worldmatPaCER(flirtmat, niiSrc, niiTrg)
if ischar(niiSrc)
    niiSrc = NiftiMod(niiSrc, 'isToBeCached', false); %HA: interpreating as filename, loading nii header 
end
if ischar(niiTrg)
    niiTrg = NiftiMod(niiTrg, 'isToBeCached', false); %HA: interpreating as filename, loading nii header  
end

% src = inv(flirtmat) * trg
% srcvox = src.mat \ inv(flirtmat) * trg.mat * trgvox
% BUT, flirt doesn't use src.mat, only absolute values of the
% scaling elements from it,
% AND, if images are not radiological, the x-axis is flipped, see:
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0810&L=FSL&P=185638
%  https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind0903&L=FSL&P=R93775
trgscl = nifti2scl(niiTrg.transformationMatrix, niiTrg.voxdim);
srcscl = nifti2scl(niiSrc.transformationMatrix, niiSrc.voxdim);
fslvoxmat = inv(srcscl) * inv(flirtmat) * trgscl;

% AND, Flirt's voxels are zero-based, while SPM's are one-based...
addone = eye(4);
%addone(:, 4) = 1; We DON'T need that it PaCER (would induce off by one
%error)
spmvoxmat = addone * fslvoxmat * inv(addone);

worldmat = niiSrc.transformationMatrix * spmvoxmat * inv(niiTrg.transformationMatrix);

%% Nested Helper
    function scl = nifti2scl(mat, voxdim)
        scl = diag([sqrt(sum(mat(1:3,1:3).^2)) 1]);
        if det(mat) > 0
            % neurological, x-axis is flipped, such that [3 2 1 0] and [0 1 2 3]
            % have the same *scaled* coordinates:
            xflip = diag([-1 1 1 1]); 
            xflip(1, 4) =voxdim(1)-1;% reflect about centre
            scl = scl * xflip;
        end
    end
end