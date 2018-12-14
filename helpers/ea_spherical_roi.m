function ea_spherical_roi(fname,mni,r,crop)
if ~exist('crop','var')
    crop=1;
end
try % try with high res bb image if subcortical ROI
    Vol=ea_load_nii([ea_space,'bb.nii']);
    Vol.img(:)=nan;
    voxmm = Vol.voxsize;
    for a=1:size(mni,1)
        X= mni(a,1); Y = mni(a,2); Z = mni(a,3);
        XYZ=[X,Y,Z,ones(length(X),1)]';
        XYZ=Vol.mat\XYZ; % to voxel space.
        XYZ=(XYZ(1:3,:)');
        
        xe = XYZ(1)-round(2*r/voxmm(1)):XYZ(1)+round(2*r/voxmm(1));
        ye = XYZ(2)-round(2*r/voxmm(2)):XYZ(2)+round(2*r/voxmm(2));
        ze = XYZ(3)-round(2*r/voxmm(3)):XYZ(3)+round(2*r/voxmm(3));
        
        
        [xx yy zz] = meshgrid(1:length(xe),1:length(ye),1:length(ze));
        S = round(sqrt((xx-2*r/voxmm(1)).^2+(yy-2*r/voxmm(2)).^2+(zz-2*r/voxmm(3)).^2)<=r/voxmm(1));
        xix=squeeze(xx(1,:,1)+round(XYZ(1)-2*r/voxmm(1)));
        yiy=squeeze(yy(:,1,1)+round(XYZ(2)-2*r/voxmm(1)));
        ziz=squeeze(zz(1,1,:)+round(XYZ(3)-2*r/voxmm(1)));
         Vol.img(xix,yiy,ziz)=S;
    end
     Vol.img(Vol.img~=1)=0;
    Vol.dt =[16,0];
    Vol.fname=fname;
    ea_write_nii(Vol);
catch % if not possible revert to whole brain t1
    
    Vol=ea_load_nii([ea_space,'t1.nii']);
    
    Vol.img(:)=nan;
    voxmm = Vol.voxsize;
    for a=1:size(mni,1)
        X= mni(a,1); Y = mni(a,2); Z = mni(a,3);
        XYZ=[X,Y,Z,ones(length(X),1)]';
        XYZ=Vol.mat\XYZ; % to voxel space.
        XYZ=(XYZ(1:3,:)');
        
        xe = XYZ(1)-round(2*r/voxmm(1)):XYZ(1)+round(2*r/voxmm(1));
        ye = XYZ(2)-round(2*r/voxmm(2)):XYZ(2)+round(2*r/voxmm(2));
        ze = XYZ(3)-round(2*r/voxmm(3)):XYZ(3)+round(2*r/voxmm(3));
        
        
        [xx yy zz] = meshgrid(1:length(xe),1:length(ye),1:length(ze));
        S = round(sqrt((xx-2*r/voxmm(1)).^2+(yy-2*r/voxmm(2)).^2+(zz-2*r/voxmm(3)).^2)<=r/voxmm(1));
        xix=squeeze(xx(1,:,1)+round(XYZ(1)-2*r/voxmm(1)));
        yiy=squeeze(yy(:,1,1)+round(XYZ(2)-2*r/voxmm(1)));
        ziz=squeeze(zz(1,1,:)+round(XYZ(3)-2*r/voxmm(1)));
        Vol.img(xix,yiy,ziz)=S;
        Vol.img(Vol.img~=1)=0;
        Vol.dt =[16,0];
        Vol.fname=fname;
        ea_write_nii(Vol);
    end
end
if crop
ea_crop_nii(fname)
end
