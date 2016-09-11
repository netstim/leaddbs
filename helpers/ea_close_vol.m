function ea_close_vol(V)

V.fname=ea_niigz(V.fname);
[~,~,ext]=fileparts(V.fname);
if strcmp(ext,'.gz')
    return
end
% closes a nifti for reading/writing (i.e. gunzipping and spm_vol
gzip(V.fname);
delete(V.fname);
