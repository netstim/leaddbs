function ea_create_surf_parcellation(parcname,erode,smooth)

if ~exist('erode', 'var')
    erode=0;
end
if ~exist('smooth', 'var')
    smooth=0;
end
labfile=ea_niigz(fullfile(ea_space,'labeling',ea_rmext(parcname)));
nii=ea_load_nii(labfile);


nii.img=round(nii.img); % make sure precision wise to map to integers
idx=unique(nii.img(:));
img_size = size(nii.img);
[X_vox, Y_vox, Z_vox] = meshgrid(1:img_size(1), 1:img_size(2), 1:img_size(3));
XYZ_mm = ea_vox2mm([X_vox(:), Y_vox(:), Z_vox(:)], nii.mat);
X = reshape(XYZ_mm(:,1), size(X_vox));
Y = reshape(XYZ_mm(:,2), size(Y_vox));
Z = reshape(XYZ_mm(:,3), size(Z_vox));

cnt=1;
ea_dispercent(0,'Iterating parcels');
for id=idx'
    if id % ignore zeroes

        onii=nii;
        onii.img=onii.img==id;
        if erode
            uid=ea_generate_uuid;
            out=fullfile(ea_getleadtempdir,[uid,'.nii']);
            onii.fname=out;
            ea_write_nii(onii);
            spm_smooth(onii.fname,fullfile(ea_getleadtempdir,['s',uid,'.nii']),[erode,erode,erode]);
            onii=ea_load_nii(fullfile(ea_getleadtempdir,['s',uid,'.nii']));
            onii.img=onii.img>0.5;
            ea_delete(onii.fname);
            ea_delete(fullfile(ea_getleadtempdir,['s',uid,'.nii']));
        end

        fv(cnt)=isosurface(X,Y,Z,permute(double(onii.img),[2,1,3]),0);
        fvc=isocaps(X,Y,Z,permute(onii.img,[2,1,3]),0);
        fv(cnt).faces=[fv(cnt).faces;fvc.faces+size(fv(cnt).vertices,1)];
        fv(cnt).vertices=[fv(cnt).vertices;fvc.vertices];
        if smooth
            fv(cnt)=ea_smoothpatch(fv(cnt),1,smooth);
        end
        cnt=cnt+1;
    end
    ea_dispercent(cnt/length(idx));
end
ea_dispercent(1,'end');
fv=ea_concatfv(fv);

plywrite(fullfile(ea_space,'labeling',[ea_rmext(parcname),'.ply']),fv.faces,fv.vertices);

ea_stlwrite(fullfile(ea_space,'labeling',[ea_rmext(parcname),'.stl']),fv);
