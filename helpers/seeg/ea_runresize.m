function status = ea_runresize(options)
    status = 0;  % default
    preproc_files = cellfun(@(x) options.subj.preproc.anat.preop.(x), ...
                            fieldnames(options.subj.preproc.anat.preop), ...
                            'uni', 0);
    MRFile = preproc_files{1};
    mmvox = 0.4;                   
    mxscale = repmat(mmvox, 1, 3);     % [0.4 0.4 0.4]

    V = spm_vol(MRFile);
    img = spm_read_vols(V);

    mrdim   = V.dim;                               % [nx ny nz]
    mrscale = sqrt(sum(V.mat(1:3,1:3).^2, 1));    % [dx dy dz] mm/vox

    mrnewdim = ceil((mrscale ./ mxscale) .* mrdim);  % target [Nx Ny Nz]
    img_res = imresize3(img, mrnewdim);

    img_info = V;

    scale_vec = mrdim ./ mrnewdim;

    Mscale = spm_matrix([0 0 0 0 0 0 scale_vec 0 0 0]);

    newmat = Mscale * img_info.mat;
    newmat(1:3,4) = img_info.mat(1:3,4);

    img_info.mat = newmat;
    img_info.dim = size(img_res);

    spm_write_vol(img_info, img_res);

    status = 1;
end