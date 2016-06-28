function mrs = DPS_nifti_to_mrs(nii)

if isstr(nii),
    nii = load_untouch_nii(nii);
end;

nii.mat = [nii.hdr.hist.srow_x ; nii.hdr.hist.srow_y ; nii.hdr.hist.srow_z; 0 0 0 1];

if ndims(nii.img) > 3
    mrs = mrstruct_init( 'series3D', double(nii.img) );
else
    mrs = mrstruct_init( 'volume', double(nii.img) );
end

mrs.vox = sqrt( sum((nii.mat*diag([1 1 1 0])).^2) );
mrs.edges = nii.mat;
mrs.edges = diag([-1 -1 1 1 ]) * mrs.edges;


%extract TR and TE
if ~isempty(nii.hdr.hk.db_name)
    tString = nii.hdr.hk.db_name;
    names = regexp(tString, '?TR:(?<TR>.+)\sTE:(?<TE>.+)', 'names');
    if ~isempty(names)
        mrs.tr = str2double(names.TR);
        mrs.te = str2double(names.TE);
    end
end
