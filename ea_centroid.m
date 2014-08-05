function stats=ea_centroid(bw)
% small replacement for the regionprops function (centroid only) which
% requires the image processing toolbox.

[xx,yy]=find(bw);
xy=[yy,xx];
stats.Centroid=mean(xy,1);