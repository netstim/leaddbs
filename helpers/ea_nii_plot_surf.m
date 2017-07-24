function gii = ea_nii_plot_surf(nii, giithresh, maxdistance)
% Plot a colored surface base on NIfTI file (heatmap, activition map...)
% The color is coded by the intensity in NIfTI file.
% 'spm_surf' will be used to extract the surface when needed.
% MATLAB GIfTI toobox is used for rendering.
%
% The second parameter 'giithresh' can either be the threshold used to
% extract the surface or the name of an already extracted GIfTI file.

if nargin < 2
    giithresh = 0.5;
end

if nargin < 3
    maxdistance = 0.5;
end

if ischar(giithresh)
    gii = giithresh;
else
    output = spm_surf(nii, 2, giithresh);
    gii = output.surffile{1};
end

surf = gifti(gii);

vol = spm_vol(nii);
img = spm_read_vols(vol);

ind = find(img);
[Xvox, Yvox, Zvox] = ind2sub(size(img), ind);
vox = [Xvox, Yvox, Zvox];
mm = single(ea_vox2mm(vox, vol.mat));

v.cdata = zeros(size(surf.vertices,1),1);
for i=1:size(surf.vertices,1)
    dist = ea_pdist2(mm, surf.vertices(i,:));
    nearest = find(dist==min(dist),1);

    if min(dist) <= maxdistance
        v.cdata(i) = img(vox(nearest,1),vox(nearest,2),vox(nearest,3));
    else
        v.cdata(i) = nan;
    end
end

figure;
plot(surf, v);
