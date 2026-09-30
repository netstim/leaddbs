function vis_cortex_seeg(~,~,resultfig,options)
% Render SEEG-style cortex and open appearance controls.


% Remove previous SEEG cortex if it already exists

oldH = getappdata(resultfig,'seegcortex');

if ~isempty(oldH)
    for s = 1:numel(oldH)
        if ~isempty(oldH{s}) && isgraphics(oldH{s})
            delete(oldH{s});
        end
    end
end

% Also clean up anything from older version
delete(findobj(resultfig,'Tag','seeg_cortex'));
delete(findobj(resultfig,'Tag','seeg_cortex_rh'));
delete(findobj(resultfig,'Tag','seeg_cortex_lh'));


% Load cortex

% cortexFile = fullfile( ...
%     ea_space(options,'cortex'), ...
%     'CortexHiRes.mat');
%
% if ~isfile(cortexFile)
%     ea_error('Missing Template Cortex');
%     return
% end
%
% S = load( ...
%     cortexFile, ...
%     'Vertices_rh','Faces_rh', ...
%     'Vertices_lh','Faces_lh');

% Load cortex

useNative = isfield(options,'native') && options.native;

if useNative

    % Patient-specific native-space cortex
    surfDir = fullfile(options.subj.subjDir, 'surfaces');

    surfFiles = dir(fullfile( ...
        surfDir, ...
        '*_space-anchorNative_*_surface.mat'));

    if isempty(surfFiles)
        ea_error(['Missing patient-specific native cortex in: ' surfDir]);
        return
    end

    % There should normally be one anchorNative surface file
    surfFile = fullfile( ...
        surfFiles(1).folder, ...
        surfFiles(1).name);

    S = load(surfFile, 'BrainSurfRaw');

    if ~isfield(S,'BrainSurfRaw') || ...
            ~isfield(S.BrainSurfRaw,'vertices') || ...
            ~isfield(S.BrainSurfRaw,'faces')

        ea_error(['Invalid patient-specific cortex file: ' surfFile]);
        return
    end

else

    % Existing MNI template cortex
    cortexFile = fullfile( ...
        ea_space(options,'cortex'), ...
        'CortexHiRes.mat');

    if ~isfile(cortexFile)
        ea_error('Missing Template Cortex');
        return
    end

    S = load( ...
        cortexFile, ...
        'Vertices_rh','Faces_rh', ...
        'Vertices_lh','Faces_lh');

end
% Default appearance

alpha = 0.4;
baseColor = [0 0 0];


% Render RIGHT and LEFT independently
% This is necessary so each hemisphere can be toggled separately.

% set(0,'CurrentFigure',resultfig);
% hold on
%
% seegH = cell(1,2);
%
% % Right hemisphere
% seegH{1} = renderHemisphere( ...
%     resultfig, ...
%     S.Vertices_rh, ...
%     S.Faces_rh, ...
%     baseColor, ...
%     alpha, ...
%     'seeg_cortex_rh');
%
% % Left hemisphere
% seegH{2} = renderHemisphere( ...
%     resultfig, ...
%     S.Vertices_lh, ...
%     S.Faces_lh, ...
%     baseColor, ...
%     alpha, ...
%     'seeg_cortex_lh');
% Render cortex

set(0,'CurrentFigure',resultfig);
hold on

seegH = cell(1,2);

if useNative

    % Patient-specific native surface is stored as one whole-brain mesh
    seegH{1} = renderHemisphere( ...
        resultfig, ...
        S.BrainSurfRaw.vertices, ...
        S.BrainSurfRaw.faces, ...
        baseColor, ...
        alpha, ...
        'seeg_cortex');

    seegH{2} = [];

else

    % MNI template: render right and left independently

    % Right hemisphere
    seegH{1} = renderHemisphere( ...
        resultfig, ...
        S.Vertices_rh, ...
        S.Faces_rh, ...
        baseColor, ...
        alpha, ...
        'seeg_cortex_rh');

    % Left hemisphere
    seegH{2} = renderHemisphere( ...
        resultfig, ...
        S.Vertices_lh, ...
        S.Faces_lh, ...
        baseColor, ...
        alpha, ...
        'seeg_cortex_lh');

end
axis vis3d
axis off

% Store handles on the main 3D figure
setappdata(resultfig,'seegcortex',seegH);


% Open the SEEG cortex controls

cwin = ea_seegcortexcontrol( ...
    resultfig, ...
    baseColor, ...
    alpha);

setappdata(resultfig,'seegcortexcontrol',cwin);

try
    WinOnTop(cwin,true);
end

end


% RENDER ONE HEMISPHERE

function p = renderHemisphere(resultfig,V,F,baseColor,alpha,tag)

% Keep full resolution
keepFraction = 1;

F = fliplr(F);

[F,V] = reducepatch(F,V,keepFraction);


% Curvature

L = cotLaplacian(V,F);

TR = triangulation(F,V);
N = vertexNormal(TR);

H = sqrt(sum((L*V).^2,2));
w_abs = abs(H);

nrmlz = @(x) min(max( ...
    (x-prctile(x,5)) ./ ...
    max(eps,prctile(x,95)-prctile(x,5)), ...
    0),1);

wc = nrmlz(w_abs).^0.3;


% Rim lighting

ax = get(resultfig,'CurrentAxes');

if isempty(ax)
    ax = gca;
end

cp = get(ax,'CameraPosition');

V2C = cp - V;
V2C = V2C ./ vecnorm(V2C,2,2);

rim = 1 - abs(sum(N .* V2C,2));
rim = rim.^3;


% Original SEEG cortex appearance

a_curve = 0;
b_rim   = 1.0;

w = min(1, a_curve*wc + b_rim*rim);

% Base color -> white rim
C = baseColor + w .* (1-baseColor);


% Patch

p = patch( ...
    'Parent',ax, ...
    'Faces',F, ...
    'Vertices',V, ...
    'FaceVertexCData',C, ...
    'FaceColor','interp', ...
    'EdgeColor','none', ...
    'FaceAlpha',alpha, ...
    'BackFaceLighting','lit', ...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.10, ...
    'DiffuseStrength',0.95, ...
    'SpecularStrength',0.03, ...
    'SpecularExponent',8, ...
    'SpecularColorReflectance',0, ...
    'Tag',tag);

% Store information needed to recolor without recomputing geometry
ud = struct();
ud.GlowWeight = w;
ud.BaseColor = baseColor;

set(p,'UserData',ud);

end



% COTANGENT LAPLACIAN

function L = cotLaplacian(V,F)

v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

e12 = v2-v1;
e23 = v3-v2;
e31 = v1-v3;

l12 = sum(e12.^2,2);
l23 = sum(e23.^2,2);
l31 = sum(e31.^2,2);

cot1 = (l12+l31-l23) ./ ...
    sqrt(max(4*l12.*l31-(l12+l31-l23).^2,eps));

cot2 = (l23+l12-l31) ./ ...
    sqrt(max(4*l23.*l12-(l23+l12-l31).^2,eps));

cot3 = (l31+l23-l12) ./ ...
    sqrt(max(4*l31.*l23-(l31+l23-l12).^2,eps));

i = [F(:,2);F(:,3);F(:,3);F(:,1);F(:,1);F(:,2)];
j = [F(:,3);F(:,2);F(:,1);F(:,3);F(:,2);F(:,1)];

s = 0.5*[cot1;cot1;cot2;cot2;cot3;cot3];

nV = size(V,1);

W = sparse(i,j,s,nV,nV);
W = (W+W.')/2;

d = -sum(W,2);

L = spdiags(d,0,nV,nV)+W;

end