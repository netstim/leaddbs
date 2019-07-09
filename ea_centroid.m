function stats=ea_centroid(bw)
% Small replacement for the regionprops function (centroid only) which
% requires the image processing toolbox.
%
% USAGE:
%
%    stats=ea_centroid(bw)
%
% INPUT:
%    bw:
%
% OUTPUT:
%    stats:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ningfei Li, Original file
%       - Daniel Duarte, Documentation

[xx,yy]=find(bw);
xy=[yy,xx];
stats.Centroid=mean(xy,1);