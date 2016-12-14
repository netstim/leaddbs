function save_mrstruct_as_nifti(mr,fname)

    nii = make_nii(mr.dataAy,mr.vox);    
    Q = diag([-1 -1 1 1]);    
    M = Q*mr.edges;
    nii.hdr.hist.srow_x = M(1,:);
    nii.hdr.hist.srow_y = M(2,:);
    nii.hdr.hist.srow_z = M(3,:);
    nii.hdr.hist.sform_code = 1;
    nii.untouch = 1;
    nii.hdr.hist.magic = 'n+1';
    save_untouch_nii(nii,fname);


