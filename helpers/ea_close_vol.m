function ea_open_vol(V)

% closes a nifti for reading/writing (i.e. gunzipping and spm_vol
gzip(V.fname);
delete(V.fname);
