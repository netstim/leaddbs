function patchTexture(patchHandle,texture)
%PATCHTEXTURE allows texture mapping onto a patch
%
% SYNOPSIS: patchTexture(patchHandle,texture)
%
% INPUT patchHandle: handle to patch object. Can be array of handles
%		texture: texture to map onto the surface defined via patch
%                Texture can be a grayscale or an RGB image
%
% OUTPUT none - the texture is mapped onto the patch identified by
%               patchHandle by reading the image at the [x,y] coordinates
%               of the vertices of the patch. The vertex coordinates are
%               rescaled so that a maximum of the image can be read.
%
% EXAMPLE
%            [xx,yy,zz] = ndgrid(-100:100,-100:100,-100:100);
%            img = xx.^2 + yy.^2 + zz.^2 < 99^2;
%            img = convn(img,ones(5,5,5)/(5*5*5));
%            p = isosurface(img,0.5);
%            ph = patch(p);
%            set(gca,'visible','off')
%            axis equal
%            texture = imread('autumn.tif');
%            patchTexture(ph,texture)
%            view(58,82)
%
% SEE ALSO patch, surface
%
% created with MATLAB ver.: 7.10.0.499 (R2010a) on Mac OS X  Version: 10.6.3 Build: 10D573
%
% created by: Jonas Dorn
% DATE: 26-May-2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TEST INPUT
if nargin < 2 || isempty(patchHandle) || isempty(texture)
    error('patchTexture needs two nonempty input arguments');
end

if ~all(ishandle(patchHandle)) || ~all(strcmp(get(patchHandle,'type'),'patch'))
    error('patchHandle needs to be the handle to a patch object')
end

if ~isnumeric(texture) || ndims(texture) > 3 || ~any(size(texture,3) == [1,3])
    error('texture needs to be an image')
end


%% COLLECT DATA

% We want an indexed image. Thus, convert to RGB if necessary, and make
% indexed image. Assume that the user knows what they're doing when they
% supply an RGB image, or an integer image
if size(texture,3) == 1
    % when converting to RGB, make sure the image is within range
    maxT = max(texture(:));
    minT = min(texture(:));
    if ~ isinteger(texture) && (maxT>1 || minT<0)
        texture = (texture-minT)/(maxT-minT);
    end
    texture = repmat(texture,[1,1,3]);
end

sizTex = size(texture);

% index
[texture,cmap] = rgb2ind(texture,256);

% loop for multiple patch handles
for ph = patchHandle(:)'
    
    % get the x- and y-coordinates of the vertices. Scale them such that they
    % cover most of the image
    vert = get(ph,'Vertices');
    % normalize xy to 0...1
    vert = vert(:,1:2);
    minV = min(vert,[],1);
    maxV = max(vert,[],1);
    vert = bsxfun(@minus,vert,minV);
    vert = bsxfun(@rdivide,vert,(maxV-minV));
    % ...and scale. Round so that vert becomes integer indices into texture
    vert = round(bsxfun(@times,vert,sizTex(1:2)-1))+1;
    
    %% COLOR THE SURFACE
    % since texture maps to the colormap, we can simply use the scaled vertex
    % x/y to read the indices we should assign to the vertices
    colorIdx = double(texture(sub2ind(sizTex(1:2),vert(:,1),vert(:,2))))+1;
    
    set(ph,'edgeColor','none','FaceColor','interp','FaceVertexCData',cmap(colorIdx,:))
    
end




