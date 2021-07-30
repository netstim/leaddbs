function ea_specsurf(varargin)

surfc = varargin{1};
color = varargin{2};
if nargin==3
    aData = varargin{3};
end

len = get(surfc,'ZData');

cd = zeros([size(len),3]);
cd(:,:,1) = color(1);

try % works if color is denoted as 1x3 array
    cd(:,:,2) = color(2);
    cd(:,:,3) = color(3);
catch % if color is denoted as gray value (1x1) only
    cd(:,:,2) = color(1);
    cd(:,:,3) = color(1);
end

cd = cd+0.01*randn(size(cd));

set(surfc,'FaceColor','interp');
set(surfc,'CData',cd);

try % for patches
    Vertices = get(surfc,'Vertices');
    cd = zeros(size(Vertices));
    cd(:) = color(1);
    set(surfc,'FaceVertexCData',cd);
end

set(surfc,'AlphaDataMapping','none');
set(surfc,'FaceLighting','phong');
set(surfc,'SpecularColorReflectance',0);
set(surfc,'SpecularExponent',10);
set(surfc,'EdgeColor','none')

if nargin==3
    set(surfc,'FaceAlpha',aData);
end
