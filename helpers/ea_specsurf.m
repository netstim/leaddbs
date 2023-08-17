function ea_specsurf(varargin)

surfc = varargin{1};

color = varargin{2};
if length(color)==1
    color = repmat(color,1,3);
end

if nargin>=3
    aData = varargin{3};
end

len = get(surfc,'ZData');

cd = zeros([size(len),3]);
cd(:,:,1) = color(1);
cd(:,:,2) = color(2);
cd(:,:,3) = color(3);

set(surfc,'AlphaDataMapping','none');

cd = cd+0.01*randn(size(cd));
set(surfc,'CData',cd);

set(surfc,'EdgeColor','none');

if nargin>=3
    set(surfc,'FaceAlpha',aData);
end

set(surfc,'FaceColor','interp');
set(surfc,'FaceLighting','phong');

if isa(surfc,'matlab.graphics.primitive.Patch') % for patches
    Vertices = get(surfc,'Vertices');
    cd = zeros(size(Vertices));
    cd = repmat(color, size(cd,1), size(cd,2)/size(color,2));
    set(surfc,'FaceVertexCData',cd);
end

if nargin>3
    switch varargin{4} % material
        case 'metal'
            set(surfc,'AmbientStrength',0.17); %0.1;
            set(surfc,'DiffuseStrength',0.4); %0.2;
            set(surfc,'SpecularColorReflectance',0);
            set(surfc,'SpecularExponent',6);
            set(surfc,'SpecularStrength',1.0);
            % met=load([ea_getearoot,'icons',filesep,'metal.mat']);
            % ea_patchtexture(surfc,met.X);
        case 'insulation'
            set(surfc,'AmbientStrength',0.17);
            set(surfc,'DiffuseStrength',0.4);
            set(surfc,'SpecularColorReflectance',1.0);
            set(surfc,'SpecularExponent',6);
            set(surfc,'SpecularStrength',0.2);
    end
else % default
    set(surfc,'AmbientStrength',0.3);
    set(surfc,'DiffuseStrength',0.4);
    set(surfc,'SpecularColorReflectance',0);
    set(surfc,'SpecularExponent',3);
    set(surfc,'SpecularStrength',0.21);
end
