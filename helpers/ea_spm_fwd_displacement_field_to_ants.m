function [] = ea_spm_fwd_displacement_field_to_ants(in_spm_field, out_ants_field)

if strcmp(out_ants_field(end-2:end), '.gz')
    out_ants_field = out_ants_field(1:end-3);
end

n = load_nii(in_spm_field);
v = spm_vol(in_spm_field);

s = n.hdr.dime.dim(2:4);
index = 1:prod(s);
[v1,v2,v3] = ind2sub(s,index);
mm = ea_vox2mm([v1',v2',v3'], v.mat);
mm = reshape(mm, [s,1,3]);
n.img = n.img - mm;

n.img(:,:,:,1,1:2) = -n.img(:,:,:,1,1:2);

n.hdr.dime.intent_code = 1007;
n.hdr.dime.xyzt_units = 2;
n.hdr.hist.descrip = '';

save_nii(n,out_ants_field,[]);
gzip(out_ants_field);
delete(out_ants_field);

end