function ea_ftr2nii(ftrFile, reference, outputName)

if ~exist('reference','var') || isempty(reference)
    reference = [ea_space,'t1.nii'];
end
refnii = ea_load_nii(reference);

if ~exist('outputName','var')
    outputName = regexprep(ftrFile, '\.mat$', '.nii');
end

[fibers,idx,voxmm,mat,vals] = ea_loadfibertracts(ftrFile);

fiberval = repelem(vals, idx);
fiberval(isnan(fiberval)) = 0;
fiberval(isinf(fiberval)) = 0;

if strcmp(voxmm, 'mm')
    % Covert fibers into the reference voxel space
    fibers_vox = round(ea_mm2vox(fibers(:,1:3), refnii.mat));
elseif strcmp(voxmm, 'vox')
    if isempty(mat) || all(mat==refnii.mat, 'all')
        % Suppose fibers are already in the reference voxel space
        fibers_vox = round(fibers(:,1:3));
    else
        % Covert fibers to mm and then into the reference voxel space
        fibers_mm = ea_vox2mm(fibers(:,1:3), mat);
        fibers_vox = round(ea_mm2vox(fibers_mm(:,1:3), refnii.mat));
    end
end

% Remove all outliers
fibers_voxc = fibers_vox(all(fibers_vox>0, 2) & all(fibers_vox<=refnii.dim, 2), :);
fibervalc = fiberval(all(fibers_vox>0, 2) & all(fibers_vox<=refnii.dim, 2), :);

% Calculate values at voxels
refnii.img = zeros(size(refnii.img));
fibImgInd = sub2ind(size(refnii.img), fibers_voxc(:,1), fibers_voxc(:,2), fibers_voxc(:,3));
[unifibImgInd, ~, ic] = unique(fibImgInd);
fibImgVal = accumarray(ic, fibervalc);
refnii.img(unifibImgInd) = refnii.img(unifibImgInd) + fibImgVal;

% Save nifti
refnii.fname = outputName;
refnii.dt(1) = 16;
ea_write_nii(refnii);
