function h=ea_parcellation2multisurface(nii,indices,colormap)

% surface can be a char (smoothed/full) or a cell pointing to rh and lh
% .mz3 surface files.


if ~exist('surface','var')
    surface='auto';
end
if ~exist('sides','var')
    sides=1:2;
end

if ~exist('colormap','var')
    colormap=lines; %ea_color_wes('all');
end
%gradientLevel = length(colormap);







if ischar(nii)
    if exist(nii,'file')
    res=ea_load_nii(nii);
    else
        try
            res=ea_load_nii(ea_niigz(fullfile(ea_space([],'labeling'),[nii,'.nii'])));
        catch
            ea_error('Input File does not exist.');
        end
    end
else
    res=nii;
end


parcels=unique(res.img(:))';

if exist('indices','var')
idx=ismember(parcels,indices);
parcels=parcels(idx);
end

while size(colormap,1)<length(parcels)
    colormap=[colormap;colormap];
end

res.img=round(res.img);
cnt=1;
img_size = size(res.img);
[X_vox, Y_vox, Z_vox] = meshgrid(1:img_size(1), 1:img_size(2), 1:img_size(3));
XYZ_mm = ea_vox2mm([X_vox(:), Y_vox(:), Z_vox(:)], res.mat);
X = reshape(XYZ_mm(:,1), size(X_vox));
Y = reshape(XYZ_mm(:,2), size(Y_vox));
Z = reshape(XYZ_mm(:,3), size(Z_vox));

for parcel=parcels
    if parcel % ignore zero
    thisresp=res.img;
    thisresp=res.img==parcel;


    pcl=isosurface(X,Y,Z,permute(thisresp,[2,1,3]),0.5);

%    pcl.vertices=res.mat*[pcl.vertices';ones(1,size(pcl.vertices,1))];
%    pcl.vertices=pcl.vertices(1:3,:)';
    pcl=ea_smoothpatch(pcl,1,3);
    try
    h{cnt}=patch('Faces',pcl.faces,'Vertices',pcl.vertices,'FaceColor',colormap(cnt,:),'EdgeColor','none',...
        'SpecularStrength',0.35,'SpecularExponent',30,'SpecularColorReflectance',0,'AmbientStrength',0.07,'DiffuseStrength',0.45,'FaceLighting','gouraud');
    catch
        keyboard
    end
    cnt=cnt+1;
    end
end

