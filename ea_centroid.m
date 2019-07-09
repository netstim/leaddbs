function stats=ea_centroid(bw)
% Small replacement for the regionprops function (centroid only) which
% requires the image processing toolbox.
%
% USAGE:
%
%    stats=ea_centroid(bw)
%
% INPUTS:
%    bw:
%
% OUTPUTS:
%    stats:
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ning Fey, Original file
%       - Daniel Duarte, Documentation

[xx,yy]=find(bw);
xy=[yy,xx];
stats.Centroid=mean(xy,1);