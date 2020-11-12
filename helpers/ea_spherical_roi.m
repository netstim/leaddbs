function ea_spherical_roi(fname,mni,r,crop,template,bg)
if ~exist('crop','var')
    crop=1;
end

if exist('template','var')
    Vol=ea_load_nii(template);
else
    Vol=ea_load_nii([ea_space,'t1.nii']);
end

if ~exist('bg','var')
    Vol.img(:)=nan;
else
    Vol.img(:)=bg;
end

voxmm = Vol.voxsize;
for a=1:size(mni,1)
    X= mni(a,1); Y = mni(a,2); Z = mni(a,3);
    XYZ=[X,Y,Z,ones(length(X),1)]';
    XYZ=Vol.mat\XYZ; % to voxel space.
    XYZ=(XYZ(1:3,:)');

    xe = XYZ(1)-round(2*r/voxmm(1)):XYZ(1)+round(2*r/voxmm(1));
    ye = XYZ(2)-round(2*r/voxmm(2)):XYZ(2)+round(2*r/voxmm(2));
    ze = XYZ(3)-round(2*r/voxmm(3)):XYZ(3)+round(2*r/voxmm(3));

    [xx, yy, zz] = meshgrid(1:length(xe),1:length(ye),1:length(ze));
    S = sqrt((xx-2*r/voxmm(1)).^2+(yy-2*r/voxmm(2)).^2+(zz-2*r/voxmm(3)).^2)<=r/voxmm(1);

    xix=squeeze(xx(1,:,1)+round(XYZ(1)-2*r/voxmm(1)))';
    yiy=squeeze(yy(:,1,1)+round(XYZ(2)-2*r/voxmm(1)));
    ziz=squeeze(zz(1,1,:)+round(XYZ(3)-2*r/voxmm(1)));

    try
        Vol.img(xix,yiy,ziz)=S;
    catch % negative indices.
        for xxx=1:length(xix)
            for yyy=1:length(yiy)
                for zzz=1:length(ziz)
                    try
                        Vol.img(xix(xxx),yiy(yyy),ziz(zzz))=S(xxx,yyy,zzz);
                    end
                end
            end
        end
    end
    Vol.img(Vol.img~=1)=0;
    Vol.dt =[16,0];
    Vol.fname=fname;
    Vol.img=Vol.img(1:Vol.dim(1),1:Vol.dim(2),1:Vol.dim(3));
    ea_write_nii(Vol);
end

if crop
    ea_crop_nii(fname)
end
