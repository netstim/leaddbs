function roi = ea_spherical_roi(fname,center,radius,crop,ref,bg)

% Write out NIfTI or not
if isempty(fname)
    writeoutNii = 0;
else
    writeoutNii = 1;
end

% Expand radius in case multiple centers specified
if size(center,1)>1
    if length(radius)==1
        radius = repmat(radius, 1, size(center,1));
    else
        error('Length of centers doesn''t match length of radius!');
    end
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

voxsize = ref.voxsize;
for i=1:size(center,1)
    % mm to voxel conversion
    XYZ = ea_mm2vox(center(i,:), ref.mat);
    r = radius(i);

    xe = XYZ(1)-round(r/voxsize(1)):XYZ(1)+round(r/voxsize(1));
    ye = XYZ(2)-round(r/voxsize(2)):XYZ(2)+round(r/voxsize(2));
    ze = XYZ(3)-round(r/voxsize(3)):XYZ(3)+round(r/voxsize(3));

    [xx, yy, zz] = meshgrid(1:length(xe),1:length(ye),1:length(ze));
    S = sqrt((xx-r/voxsize(1)).^2+(yy-r/voxsize(2)).^2+(zz-r/voxsize(3)).^2)<=r/mean(voxsize);

    xix=squeeze(xx(1,:,1)+round(XYZ(1)-r/voxsize(1)-1))';
    yiy=squeeze(yy(:,1,1)+round(XYZ(2)-r/voxsize(2)-1));
    ziz=squeeze(zz(1,1,:)+round(XYZ(3)-r/voxsize(3)-1));

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
        ea_autocrop(fname)
    end
end
