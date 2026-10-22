function TransformFixedParameters = ea_field_ref2itk(reference_image)
%EA_FIELD_REF2ITK Build ITK displacement-field fixed parameters from a NIfTI grid.

hdr = ea_fslhd(reference_image);

spacing = [hdr.pixdim1, hdr.pixdim2, hdr.pixdim3];
if any(~isfinite(spacing)) || any(spacing <= 0)
    error('Invalid voxel spacing in reference image: %s', reference_image);
end

stoX = hdr.sto_xyz1(1:3);
stoY = hdr.sto_xyz2(1:3);
stoZ = hdr.sto_xyz3(1:3);
voxToRAS = [stoX(:)'; stoY(:)'; stoZ(:)'];

rasToLPS = diag([-1, -1, 1]);
direction = rasToLPS * voxToRAS * diag(1 ./ spacing);

TransformFixedParameters = zeros(18, 1);
TransformFixedParameters(1:3) = [hdr.dim1; hdr.dim2; hdr.dim3];
TransformFixedParameters(4:6) = [-hdr.sto_xyz1(4); -hdr.sto_xyz2(4); hdr.sto_xyz3(4)];
TransformFixedParameters(7:9) = spacing(:);
TransformFixedParameters(10:18) = reshape(direction', [], 1);

end
