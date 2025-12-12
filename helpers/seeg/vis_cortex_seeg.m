try delete(cortexH{1}); end
try delete(cortexH{2}); end
delete(findobj(gca, 'Tag', 'seeg_cortex'));

S = load('CortexHiRes.mat','Vertices_rh','Faces_rh','Vertices_lh','Faces_lh');

Vfull = [S.Vertices_lh; S.Vertices_rh];
Ffull = [fliplr(S.Faces_lh); fliplr(S.Faces_rh) + size(S.Vertices_lh,1)];
% If you only want one hemisphere
% Vfull = [S.Vertices_lh];
% Ffull = [fliplr(S.Faces_lh)];

keepFraction = 1;                  % resolution, 0-1

[F, V] = reducepatch(Ffull, Vfull, keepFraction);

TR = triangulation(F,V);
A  = TR.freeBoundary;

% Geometric transforms for 'inverse' cortex visualization
L  = cotLaplacian(V,F);  
TR = triangulation(F,V);
N  = vertexNormal(TR);    
H     = sqrt(sum((L * V).^2, 2));
w_abs = abs(H);
nrmlz = @(x) min(max((x - prctile(x,5)) ./ max(eps, prctile(x,95)-prctile(x,5)),0),1);
wc    = nrmlz(w_abs).^.3;

% Lighting manipulations for 'darker' faces and 'glow' around the rim
ax  = gca;
cp  = get(ax,'CameraPosition');
V2C = cp - V;  V2C = V2C ./ vecnorm(V2C,2,2);
rim = 1 - abs(sum(N .* V2C, 2));
rim = rim.^2;

figure(gcf); hold on;

% Parameters for cortex glow
a_curve = 0.1;   % glow of the face / outer parts of the cortex, 0-1
b_rim   = 1.0;   % glow of the rim, 0-1
w = min(1, a_curve*wc + b_rim*rim);

% Dark base + glow to white
base = [0 0 0];
C = base + w .* (1 - base);

% Render
p = patch('Faces',F,'Vertices',V, ...
          'FaceVertexCData',C, 'FaceColor','interp', ...
          'EdgeColor','none', 'FaceAlpha',0.4, ...   % slightly more transparent
          'BackFaceLighting','lit', 'FaceLighting','gouraud', ...
          'AmbientStrength',0.10, 'DiffuseStrength',0.95, ...
          'SpecularStrength',0.03, 'SpecularExponent',8, ...
          'SpecularColorReflectance',0, 'Tag', 'seeg_cortex');
axis vis3d off

% Math helpers
function L = cotLaplacian(V,F)
  v1 = V(F(:,1),:); v2 = V(F(:,2),:); v3 = V(F(:,3),:);
  e12 = v2 - v1; e23 = v3 - v2; e31 = v1 - v3;

  l12 = sum(e12.^2,2); l23 = sum(e23.^2,2); l31 = sum(e31.^2,2);

  cot1 = (l12 + l31 - l23) ./ sqrt(max(4*l12.*l31 - (l12 + l31 - l23).^2, eps));
  cot2 = (l23 + l12 - l31) ./ sqrt(max(4*l23.*l12 - (l23 + l12 - l31).^2, eps));
  cot3 = (l31 + l23 - l12) ./ sqrt(max(4*l31.*l23 - (l31 + l23 - l12).^2, eps));

  i = [F(:,2); F(:,3); F(:,3); F(:,1); F(:,1); F(:,2)];
  j = [F(:,3); F(:,2); F(:,1); F(:,3); F(:,2); F(:,1)];
  s = 0.5 * [cot1; cot1; cot2; cot2; cot3; cot3];
  nV = size(V,1);
  W = sparse(i, j, s, nV, nV); W = (W + W.')/2;

  d = -sum(W,2);
  L = spdiags(d,0,nV,nV) + W;
end
