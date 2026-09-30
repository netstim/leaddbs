function [hIns, hCon, projPts] = ea_render_curved_lead(Ct, coords, elspec, opt)
%EA_RENDER_CURVED_LEAD  Curved DBS lead rendered as cylinders:
%   - Contacts: cylinders centered on projected coords
%   - Insulation: cylinders ONLY in the gaps between contacts (no overlap)
%
% This avoids z-fighting entirely and makes insulation render like a solid surface.

% ---------- defaults ----------
if nargin < 3 || isempty(elspec), elspec = struct; end
if nargin < 4 || isempty(opt),    opt    = struct; end

if ~isfield(elspec,'lead_diameter'),    elspec.lead_diameter    = 1.27; end
if ~isfield(elspec,'contact_diameter'), elspec.contact_diameter = elspec.lead_diameter; end
if ~isfield(elspec,'contact_height'),   elspec.contact_height   = 1.5;  end

% --- LOOK / STYLE ---
% Match ea_showelectrode's DBS material convention: insulation uses the
% electrode lead color and contacts use the metallic contact color.
if ~isfield(elspec,'insulation_alpha'), elspec.insulation_alpha = 1.0; end
if ~isfield(elspec,'lead_color'),       elspec.lead_color       = [0.8 0.8 0.8]; end
if ~isfield(elspec,'insulation_color'), elspec.insulation_color = elspec.lead_color; end
if ~isfield(elspec,'contact_color'),    elspec.contact_color    = [0.7 0.7 0.7]; end
if ~isfield(opt,'active_idx'),      opt.active_idx = []; end
if ~isfield(opt,'active_color'),    opt.active_color = [0.8 0.2 0.2]; end
if ~isfield(opt,'contact_alpha'),   opt.contact_alpha = 1; end

if ~isfield(opt,'nTheta'),          opt.nTheta = 24; end     % insulation cylinder resolution
if ~isfield(opt,'contact_nTheta'),  opt.contact_nTheta = 18; end

if ~isfield(opt,'subsample_step'),  opt.subsample_step = 1; end

% geometry robustness
if ~isfield(opt,'min_centerline_step'), opt.min_centerline_step = 0.05; end % mm
if ~isfield(opt,'resample_step'),       opt.resample_step       = 0.25; end % mm

% lighting
if ~isfield(opt,'faceLighting'),     opt.faceLighting = 'phong'; end
if ~isfield(opt,'specularStrength'), opt.specularStrength = 0.2; end
if ~isfield(opt,'contactSpecularStrength'), opt.contactSpecularStrength = 1.0; end

% optional: tiny "gap" so insulation doesn't touch contact endcaps (prevents boundary flicker)
if ~isfield(opt,'end_gap_mm'), opt.end_gap_mm = 0.02; end  % mm

% --- sanitize colors ---
elspec.contact_color    = sanitize_rgb(elspec.contact_color,    [0.70 0.70 0.70]);
elspec.insulation_color = sanitize_rgb(elspec.insulation_color, [0.80 0.80 0.80]);
opt.active_color        = sanitize_rgb(opt.active_color,        [0.80 0.20 0.20]);

Ct     = double(Ct);
coords = double(coords);

% Optional subsample (coarse)
Ct = Ct(1:opt.subsample_step:end, :);
if size(Ct,1) < 2, error('Ct too short'); end

% ---------- clean + resample Ct ----------
Ct = remove_tiny_steps(Ct, opt.min_centerline_step);
if size(Ct,1) < 3, error('Ct too short after cleaning'); end

Ct = resample_centerline_arclength(Ct, opt.resample_step);
Ct = remove_tiny_steps(Ct, opt.min_centerline_step);
if size(Ct,1) < 4, error('Ct too short after resampling'); end

% ---------- radii/heights ----------
r_ins  = elspec.lead_diameter/2;
r_cont = (elspec.contact_diameter/2);  % can be same as shaft since we no longer overlap
h_cont = elspec.contact_height;


% ---------- distal tip ----------
if isfield(elspec,'tip_length') && ~isempty(elspec.tip_length)
    tip_length = double(elspec.tip_length);
else
    tip_length = 0;
end

% ---------- RMF frame along Ct (stable tangents) ----------
[T, N, B] = rmf_frame(Ct); %#ok<ASGLU>  % N,B not required here but kept if you need later

% ---------- arc-length parameterization of Ct ----------
[sCt, segLen] = centerline_arclength(Ct);
sEnd = sCt(end);

% Start cap
[P_start, t_start] = eval_centerline_at_s(Ct, sCt, 0);

% End cap
[P_end,   t_end]   = eval_centerline_at_s(Ct, sCt, sEnd);

% ---------- project coords -> polyline (get arclength positions) ----------
M = size(coords,1);
projPts = zeros(M,3);
sM = zeros(M,1);      % arclength position of each projected coordinate

if M > 0
    A  = Ct(1:end-1,:);
    Bp = Ct(2:end,:);
    AB = Bp - A;
    AB2 = sum(AB.^2,2) + eps;

    for m = 1:M
        P = coords(m,:);
        [Pm, segIdx, tseg] = project_point_to_polyline(P, A, AB, AB2);
        projPts(m,:) = Pm;
        sM(m) = sCt(segIdx) + tseg * segLen(segIdx);
    end
else
    hIns = gobjects(0,1);
    hCon = gobjects(0,1);
    return;
end

% Sort contacts along lead by arclength (important for “between contacts” gaps)
[sM_sorted, ord] = sort(sM, 'ascend');
projPts_sorted = projPts(ord,:);

% Preallocate handles
hCon_sorted = gobjects(M,1);

% For each contact: compute start/end along CURVE using arclength, not straight tangent
halfH = 0.5*h_cont;
gap   = opt.end_gap_mm;

sStart = max(0,     sM_sorted - halfH);
sStop  = min(sEnd,  sM_sorted + halfH);

% Also shrink start/stop slightly to leave a tiny gap to insulation segments (optional)
sStart = min(sStop, sStart + gap);
sStop  = max(sStart, sStop - gap);

P0 = zeros(M,3);
P1 = zeros(M,3);

%distal electrode tip
if tip_length > 0

    % Distal edge of first contact along centerline
    sTipBase = sStart(1);

    % Direction pointing outward from the first contact
    [P_tipBase, t_tip] = eval_centerline_at_s(Ct, sCt, sTipBase);

    % Point tip direction AWAY from the rest of the electrode
    t_tip = -t_tip;

    % Build rounded tip extending tip_length beyond first contact
    hTip = local_rounded_tip( ...
        P_tipBase, ...
        t_tip, ...
        tip_length, ...
        r_ins, ...
        opt.nTheta);

    set(hTip, ...
        'FaceColor', elspec.insulation_color, ...
        'FaceAlpha', elspec.insulation_alpha, ...
        'EdgeColor', 'none', ...
        'FaceLighting', opt.faceLighting, ...
        'SpecularColorReflectance', 1.0, ...
        'SpecularExponent', 6, ...
        'SpecularStrength', opt.specularStrength, ...
        'AmbientStrength', 0.17, ...
        'DiffuseStrength', 0.4, ...
        'Tag', 'LeadTip');
end

for k = 1:M
    [P0(k,:), ~] = eval_centerline_at_s(Ct, sCt, sStart(k));
    [P1(k,:), ~] = eval_centerline_at_s(Ct, sCt, sStop(k));

    [Pc, ~] = eval_centerline_at_s(Ct, sCt, sM_sorted(k));
    projPts_sorted(k,:) = Pc;                 % ensure contact center is on Ct exactly

    % Build contact cylinder by endpoints on the curve
    hc = local_oriented_cylinder(P0(k,:), P1(k,:), r_cont, opt.contact_nTheta);

    isActive = ismember(ord(k), opt.active_idx);
    set(hc, 'FaceColor', tern(isActive, opt.active_color, elspec.contact_color), ...
            'FaceAlpha', opt.contact_alpha, ...
            'EdgeColor', 'none', ...
            'FaceLighting', opt.faceLighting, ...
            'SpecularColorReflectance', 0, ...
            'SpecularExponent', 6, ...
            'SpecularStrength', opt.contactSpecularStrength, ...
            'AmbientStrength', 0.17, ...
            'DiffuseStrength', 0.4, ...
            'Tag', sprintf('Contact_%d', ord(k)));

    hCon_sorted(k) = hc;
end

% Unsort outputs back to original order
hCon = gobjects(M,1);
hCon(ord) = hCon_sorted;
projPts = zeros(M,3);
projPts(ord,:) = projPts_sorted;

% ---------- insulation segments BETWEEN contacts ----------
% We create cylinders on the curve for:
%   [0 -> first contact start], [end of contact i -> start of contact i+1], [last contact end -> sEnd]
segS0 = [0;      sStop];   % (M+1)x1
segS1 = [sStart; sEnd];    % (M+1)x1

% remove tiny/negative segments (also protects against overlapping contacts)
minLen = (2*gap + 1e-3);    % mm
keep = (segS1 - segS0) > minLen;

segS0 = segS0(keep);
segS1 = segS1(keep);

hIns = gobjects(numel(segS0),1);

for j = 1:numel(segS0)
    [Q0, ~] = eval_centerline_at_s(Ct, sCt, segS0(j) + gap);
    [Q1, ~] = eval_centerline_at_s(Ct, sCt, segS1(j) - gap);

    if norm(Q1 - Q0) < 1e-3
        continue;
    end

    hs = local_oriented_cylinder(Q0, Q1, r_ins, opt.nTheta);

    set(hs, 'FaceColor', elspec.insulation_color, ...
            'FaceAlpha', elspec.insulation_alpha, ...
            'EdgeColor', 'none', ...
            'FaceLighting', opt.faceLighting, ...
            'SpecularColorReflectance', 1.0, ...
            'SpecularExponent', 6, ...
            'SpecularStrength', opt.specularStrength, ...
            'AmbientStrength', 0.17, ...
            'DiffuseStrength', 0.4, ...
            'Tag', sprintf('InsulationSeg_%d', j));

    hIns(j) = hs;
end

% ---------- end caps ----------
hCapStart = local_disk_cap(P_start, -t_start, r_ins, opt.nTheta);
set(hCapStart, ...
    'FaceColor', elspec.insulation_color, ...
    'FaceAlpha', elspec.insulation_alpha, ...
    'FaceLighting', opt.faceLighting, ...
    'SpecularColorReflectance', 1.0, ...
    'SpecularExponent', 6, ...
    'SpecularStrength', opt.specularStrength, ...
    'AmbientStrength', 0.17, ...
    'DiffuseStrength', 0.4, ...
    'Tag', 'LeadCap_Start');

hCapEnd = local_disk_cap(P_end, t_end, r_ins, opt.nTheta);
set(hCapEnd, ...
    'FaceColor', elspec.insulation_color, ...
    'FaceAlpha', elspec.insulation_alpha, ...
    'FaceLighting', opt.faceLighting, ...
    'SpecularColorReflectance', 1.0, ...
    'SpecularExponent', 6, ...
    'SpecularStrength', opt.specularStrength, ...
    'AmbientStrength', 0.17, ...
    'DiffuseStrength', 0.4, ...
    'Tag', 'LeadCap_End');

end


% =========================
% Helpers
% =========================

function Ct = remove_tiny_steps(Ct, minStep)
if size(Ct,1) < 2, return; end
d = vecnorm(diff(Ct,1,1),2,2);
keep = [true; d > minStep];
Ct = Ct(keep,:);
end

function C = resample_centerline_arclength(Ct, step)
s = [0; cumsum(vecnorm(diff(Ct,1,1),2,2))];
sEnd = s(end);
if sEnd < step || size(Ct,1) < 3
    C = Ct; return;
end
sq = (0:step:sEnd).';
C = interp1(s, Ct, sq, 'pchip');
end

function [T, N, B] = rmf_frame(Ct)
Npts = size(Ct,1);

Tp = zeros(Npts,3);
Tp(1,:)        = Ct(2,:) - Ct(1,:);
Tp(end,:)      = Ct(end,:) - Ct(end-1,:);
Tp(2:end-1,:)  = Ct(3:end,:) - Ct(1:end-2,:);
Tp = Tp ./ (vecnorm(Tp,2,2) + eps);
T = Tp;

t0 = T(1,:).';
cand = [1;0;0];
if abs(dot(cand,t0)) > 0.9, cand = [0;1;0]; end
n0 = cand - t0*dot(cand,t0);
n0 = n0/(norm(n0)+eps);

N = zeros(Npts,3);
B = zeros(Npts,3);
N(1,:) = n0.';
B(1,:) = cross(T(1,:), N(1,:));
B(1,:) = B(1,:)/(norm(B(1,:))+eps);

for i = 2:Npts
    t_prev = T(i-1,:).';
    t_curr = T(i,:).';

    ax = cross(t_prev, t_curr);
    s  = norm(ax);
    c  = dot(t_prev, t_curr);

    Ni = N(i-1,:).';
    if s > 1e-10
        ax = ax / s;
        ang = atan2(s, c);
        Rpt = local_axis_angle_to_rotm(ax, ang);
        Ni = Rpt * Ni;
    end

    Ni = Ni - t_curr*dot(t_curr,Ni);
    Ni = Ni/(norm(Ni)+eps);

    Bi = cross(t_curr, Ni);
    Bi = Bi/(norm(Bi)+eps);

    Ni = cross(Bi, t_curr);
    Ni = Ni/(norm(Ni)+eps);

    N(i,:) = Ni.';
    B(i,:) = Bi.';
end
end

function [sCt, segLen] = centerline_arclength(Ct)
segLen = vecnorm(diff(Ct,1,1),2,2);
sCt = [0; cumsum(segLen)];
end

function [P, tHat] = eval_centerline_at_s(Ct, sCt, sQuery)
% Evaluate position + tangent on polyline Ct at arclength sQuery
sQuery = max(0, min(sCt(end), sQuery));

% find segment index i such that sCt(i) <= sQuery <= sCt(i+1)
i = find(sCt <= sQuery, 1, 'last');
if i >= numel(sCt), i = numel(sCt)-1; end
if i < 1, i = 1; end

s0 = sCt(i);
s1 = sCt(i+1);
t = (sQuery - s0) / (s1 - s0 + eps);

A = Ct(i,:);
B = Ct(i+1,:);
P = A + t*(B - A);

tHat = (B - A);
tHat = tHat / (norm(tHat)+eps);
end

function [Pm, segIdx, tBest] = project_point_to_polyline(P, A, AB, AB2)
PA = A - P;
t  = -sum(PA .* AB, 2) ./ AB2;
t  = max(0, min(1, t));
Q  = A + AB .* t;
d2 = sum((Q - P).^2, 2);
[~, segIdx] = min(d2);
Pm = Q(segIdx,:);
tBest = t(segIdx);
end

function h = local_oriented_cylinder(P0, P1, r, nTheta)
[XC, YC, ZC] = cylinder(r, nTheta);

vraw = (P1 - P0);
len = norm(vraw);
if len < 1e-6
    len = 1e-6;
    v = [0 0 1];
else
    v = vraw / len;
end

ZC = ZC * len;

R = local_rotmat_from_a_to_b([0 0 1], v);

pts = [XC(:), YC(:), ZC(:)] * R.';
pts = pts + P0;

X = reshape(pts(:,1), size(XC));
Y = reshape(pts(:,2), size(YC));
Z = reshape(pts(:,3), size(ZC));
h = surf(X, Y, Z);
end

function R = local_rotmat_from_a_to_b(a, b)
a = a(:)/norm(a); b = b(:)/norm(b);
v = cross(a,b); c = dot(a,b);

if c > 1-1e-12
    R = eye(3); return
end

if c < -1+1e-12
    tmp = [1;0;0]; if abs(dot(tmp,a))>0.9, tmp=[0;1;0]; end
    axis = cross(a,tmp); axis = axis/norm(axis);
    R = local_axis_angle_to_rotm(axis, pi); return
end

s = norm(v);
vx = [   0  -v(3)  v(2);
      v(3)    0   -v(1);
     -v(2)  v(1)    0  ];
R = eye(3) + vx + vx*vx*((1-c)/(s^2 + eps));
end

function R = local_axis_angle_to_rotm(axis, th)
axis = axis(:)/norm(axis + eps);
x=axis(1); y=axis(2); z=axis(3); c=cos(th); s=sin(th); C=1-c;
R=[x*x*C+c, x*y*C-z*s, x*z*C+y*s; ...
   y*x*C+z*s, y*y*C+c, y*z*C-x*s; ...
   z*x*C-y*s, z*y*C+x*s, z*z*C+c];
end

function out = tern(cond, a, b)
if cond, out=a; else, out=b; end
end

function c = sanitize_rgb(c, fallback)
if isempty(c), c = fallback; return; end
if ischar(c) || isstring(c), c = fallback; return; end
if isnumeric(c) && size(c,2)==3 && size(c,1)>1, c = c(1,:); end
if ~isnumeric(c) || numel(c)~=3 || any(~isfinite(c)) || any(c<0) || any(c>1)
    c = fallback;
end
c = double(c(:).');
end

function h = local_disk_cap(Pc, tHat, r, nTheta)
% Draw a circular disk centered at Pc, normal to tHat
% Pc: 1x3, tHat: 1x3 or 3x1

Pc   = double(Pc(:).');                 % 1x3
tHat = double(tHat(:)); 
tHat = tHat / (norm(tHat) + eps);       % 3x1

% Orthonormal basis (u,v) spanning plane orthogonal to tHat
up = [0;0;1];
if abs(dot(up, tHat)) > 0.9
    up = [0;1;0];
end

u = cross(tHat, up);  u = u / (norm(u) + eps);
v = cross(tHat, u);   v = v / (norm(v) + eps);

theta = linspace(0, 2*pi, nTheta+1);
theta(end) = [];

ct = cos(theta(:));                     % nTheta x 1
st = sin(theta(:));                     % nTheta x 1

% Build circle points: nTheta x 3
circle = (ct*u.' + st*v.') * r;         % (n x 3)
circle = circle + repmat(Pc, nTheta, 1);

% Vertices: center first, then circle
V = [Pc; circle];                       % (nTheta+1) x 3

% Faces: triangle fan around center vertex 1
i2 = (2:(nTheta+1)).';
i3 = [3:(nTheta+1) 2].';
F  = [ones(nTheta,1) i2 i3];

h = patch('Vertices', V, 'Faces', F, 'EdgeColor', 'none');
end

function h = local_rounded_tip(Pbase, dirHat, tipLength, r, nTheta)
%LOCAL_ROUNDED_TIP
% Creates an electrode tip extending outward from Pbase.

Pbase  = double(Pbase(:).');
dirHat = double(dirHat(:).');
dirHat = dirHat / (norm(dirHat) + eps);

% Radius cannot consume more than the entire requested tip
roundLen = min(r, tipLength);

% Straight portion before rounded end
cylLength = max(0, tipLength - roundLen);

% ---- cylindrical portion ----

handles = gobjects(0);

if cylLength > 1e-6

    PcylEnd = Pbase + dirHat * cylLength;

    hcyl = local_oriented_cylinder( ...
        Pbase, ...
        PcylEnd, ...
        r, ...
        nTheta);

    handles(end+1) = hcyl;

else
    PcylEnd = Pbase;
end

% ---- rounded distal end ----
%
% Construct hemisphere whose base has radius r and whose pole
% lies exactly tipLength away from Pbase.

nPhi = max(8, round(nTheta/2));

theta = linspace(0, 2*pi, nTheta+1);
phi   = linspace(0, pi/2, nPhi);

[TH, PH] = meshgrid(theta, phi);

% Local hemisphere:
% z = 0 at base
% z = roundLen at rounded distal pole
X = r .* cos(TH) .* cos(PH);
Y = r .* sin(TH) .* cos(PH);
Z = roundLen .* sin(PH);

% Rotate local +Z direction onto electrode tip direction
R = local_rotmat_from_a_to_b([0 0 1], dirHat);

pts = [X(:), Y(:), Z(:)] * R.';
pts = pts + PcylEnd;

Xw = reshape(pts(:,1), size(X));
Yw = reshape(pts(:,2), size(Y));
Zw = reshape(pts(:,3), size(Z));

hround = surf(Xw, Yw, Zw);

handles(end+1) = hround;

% Return all tip surfaces together
h = handles;

end