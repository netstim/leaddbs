function roi = ea_spherical_roi(fname,center,radius,crop,ref,bg)

if isempty(fname)
    writeoutNii = 0;
else
    writeoutNii = 1;
end
if size(center,1)>1 && length(radius)==1
    radius = repmat(radius, 1, size(center,1));
end

% Crop the generate ROI image or not
if ~exist('crop','var')
    crop=1;
end

% Reference template image, use MNI t1 by default
if exist('template','var')
    ref = ea_load_nii(ref);
else
    ref = ea_load_nii([ea_space,'t1.nii']);
end

% Preset background
if ~exist('bg','var')
    ref.img(:)=nan;
else
    ref.img(:)=bg;
end

voxmm = ref.voxsize;
for i=1:size(center,1)
    % mm to voxel conversion
    XYZ = ea_mm2vox(center(i,:), ref.mat);
    r = radius(i);

    xe = XYZ(1)-round(2*r/voxmm(1)):XYZ(1)+round(2*r/voxmm(1));
    ye = XYZ(2)-round(2*r/voxmm(2)):XYZ(2)+round(2*r/voxmm(2));
    ze = XYZ(3)-round(2*r/voxmm(3)):XYZ(3)+round(2*r/voxmm(3));

    [xx, yy, zz] = meshgrid(1:length(xe),1:length(ye),1:length(ze));
    S = sqrt((xx-2*r/voxmm(1)).^2+(yy-2*r/voxmm(2)).^2+(zz-2*r/voxmm(3)).^2)<=r/voxmm(1);

    xix=squeeze(xx(1,:,1)+round(XYZ(1)-2*r/voxmm(1)))';
    yiy=squeeze(yy(:,1,1)+round(XYZ(2)-2*r/voxmm(1)));
    ziz=squeeze(zz(1,1,:)+round(XYZ(3)-2*r/voxmm(1)));

    try
        ref.img(xix,yiy,ziz)=S;
    catch % negative indices.
        for xxx=1:length(xix)
            for yyy=1:length(yiy)
                for zzz=1:length(ziz)
                    try
                        ref.img(xix(xxx),yiy(yyy),ziz(zzz))=S(xxx,yyy,zzz);
                    end
                end
            end
        end
    end
end

% Set ROI NIfTI structure
ref.img(ref.img~=1) = 0;
ref.dt = [16,0];
ref.img = ref.img(1:ref.dim(1),1:ref.dim(2),1:ref.dim(3));
ref.fname = fname;

roi = ref;

% Write out NIfTI
if writeoutNii
    ea_write_nii(ref);
    % Crop ROI image
    if crop
        ea_crop_nii(fname)
    end
end
