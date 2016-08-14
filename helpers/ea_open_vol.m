function V=ea_open_vol(fname)

% opens a nifti for reading/writing (i.e. gunzipping and spm_vol

if strcmp(fname(end-2:end),'.gz')
    gunzip(fname);
    delete(fname);
    fname=fname(1:end-3);
end
V=spm_vol(fname);
