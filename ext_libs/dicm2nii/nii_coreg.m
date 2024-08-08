function [M, mss] = nii_coreg(niiT, niiM, varargin)
% Affine registration of images with the same modality.
% 
% Examples:
%   M = nii_coreg(niiT, niiM)
%   M = nii_coreg('mni_T1', niiM)
%   M = nii_coreg(sub1_ses1, sub1_ses2, 'FWHM', [6 6], 'DoF', 6, 'disp', 2)
% 
% Input (first two are required): 
%   niiT: template nii name or struct, or 'mni_T1', 'mni_T2', 'mni_EPI' 
%   niiM: moving nii name or struct, same modality as the template
%   opts: options in a struct or in key/value pairs:
%     DoF: degree of freedom, default 12. Possible values:
%        3: translation (mm) only in xyz order
%        6: plus 3 rotation: rigid body alignment for the same brain
%        7: plus single scaling in xyz
%        9: from DoF=6, plus independent scaling in xyz
%       12: shearing also allowed for typical inter-subject registration
%     FWHM: FWHM of Gaussian smoothing for niiT and niiM, default [0 6] mm
%     mask: if provided, perform masked alignment after whole image registration:
%        '': no masked alignment, default unless niiT is one of the 'mni_*' 
%        'mni': use MNI brain for masked alignment. Default for 'mni_*'
%        nii: nii as mask for template (same sform/qform as niiT)
%     disp: control the display of result or progress
%        0: default, suppressing command output and viewer. 
%        1: display some information in Command Window
%        2: also show final registration result in nii_viewer
%        3: also show progress of registration for each iteration
%  (The following options are rarely needed:)
%     init: initial M (4x4) or a vector of length 3,6,7,9 or 12, matching DoF.
%        In case of poor alignment, an init guess will help
%     delta: the stop criterion used for iterations, default 1e-4.
%        It may help to reduce it if the result is very off.
%     res: resolution to take samples from niiT, default 3 mm
% 
% Output:
%   M: transformation matrix, so that
%       R = M * RM; % will transform niiM to niiT, and 
%       R = M \ RT; % will transform niiT to niiM
%     For example:
%       M = nii_coreg(niiT, niiM);
%       R = M * [niiM.hdr.sform_mat; 0 0 0 1];
%       niiM.hdr.sform_mat = R(1:3,:); % update niiM xform matrix
%       nii_viewer(niiT, niiM); % overlay to check goodness of registration
%   mss: goodness of the alignment (<0.3 typically good for mni masked alignment)

% 2305xx Xiangrui.Li@gmail.com start it from spm_affreg.m

% deal with different input argument
opts = struct('DoF', 12, 'FWHM', [0 6], 'mask', '', 'init', [], 'disp', 0, ...
    'delta', 1e-4, 'res', 3);
tkDir = fileparts(which('nii_viewer.m'));
if isequal(niiT, 'mni_T1')
    niiT = [tkDir '/templates/MNI_2mm_T1.nii'];
    opts.mask = 'mni';
elseif isequal(niiT, 'mni_T2')
    niiT = [tkDir '/templates/MNI_2mm_T2.nii'];
    opts.mask = 'mni';
elseif isequal(niiT, 'mni_EPI')
    load([tkDir '/example_data.mat'], 'MNI_2mm_EPI');
    niiT = MNI_2mm_EPI;
    opts.mask = 'mni';
end
if nargin > 2
    if nargin == 3 % struct input
        in3 = varargin{1};
        if ~isstruct(in3), error('struct or key/value pair expected'); end
    elseif nargin > 3 % key/value pairs
        in3 = struct(varargin{:});
    end
    flds = fieldnames(in3);
    for i = 1:numel(flds)
        if isfield(opts, flds{i}), opts.(flds{i}) = in3.(flds{i}); end
    end
end

if ~isstruct(niiT), niiT = nii_tool('load', niiT); end
if ~isstruct(niiM), niiM = nii_tool('load', niiM); end
if any(niiT.hdr.dim(5:end)>1), niiT.img = niiT.img(:,:,:,1); end
if any(niiM.hdr.dim(5:end)>1), niiM.img = niiM.img(:,:,:,1); end

res = max(opts.FWHM(2)/2, opts.res);
slope = niiT.hdr.scl_slope; if slope==0, slope = 1; end
imgT = double(niiT.img) * slope + niiT.hdr.scl_inter;
if opts.FWHM(1)>0
    sz = round(opts.FWHM(1) / max(niiT.hdr.pixdim(2:4)))*2 + 1;
    imgT = smooth3(imgT, 'gaussian', sz, opts.FWHM(1)/2.355);
end
[vT, dT, XT, RT, p0] = sampleT(niiT.hdr, imgT, res, opts.DoF);
muT = mean(vT);

if isequal(opts.mask, 'mni')
    load([tkDir '/example_data.mat'], 'MNI_2mm_mask');
    opts.mask = MNI_2mm_mask;
end
if ~isempty(opts.mask)
    if any(abs(niiT.hdr.sform_mat(:) - opts.mask.hdr.sform_mat(:)) > 1e-4)
        opts.mask = nii_xform(opts.mask, niiT);
    end
end

slope = niiM.hdr.scl_slope; if slope==0, slope = 1; end
imgM = double(niiM.img) * slope + niiM.hdr.scl_inter;
if opts.FWHM(2)>0
    sz = round(opts.FWHM(2) / max(niiM.hdr.pixdim(2:4)))*2 + 1;
    imgM = smooth3(imgM, 'gaussian', sz, opts.FWHM(2)/2.355);
end
F = griddedInterpolant(imgM, 'linear', 'none');

niiM.hdr.sform_code = niiT.hdr.sform_code;
RM = [niiM.hdr.sform_mat; 0 0 0 1] * [eye(4,3) [-1 -1 -1 1]'];
M = opts.init;
if isempty(M)
    cog = nii_viewer('func_handle', 'img_cog');
    cT = RT*[cog(imgT); 1] - RM*[cog(imgM); 1]; % center diff
    M = [eye(3) cT(1:3); 0 0 0 1];
elseif isvector(M)
    M = affine_mat(M(1:min(numel(M),opts.DoF)));
elseif isequal(size(M), [3 4])
    M = [M; 0 0 0 1];
elseif ~isequal(size(M), [4 4])
    error('Invalid init value for M.');
end

dB2 = opts.delta * opts.DoF;
for iter = 1:256
    IM = M * RM \ XT;
    vM = F(IM(1,:), IM(2,:), IM(3,:));
    ind = ~isnan(vM);
    a = [vM(ind); dT(:,ind)];
    b = (a * a') \ (a * vT(ind)');
    M = affine_mat(p0+b(2:end)) * M;

    if b(2:end)'*b(2:end) < dB2 % stop if little change
        if isempty(opts.mask), break; end % masked align done or not needed
        imgT(~logical(opts.mask.img)) = 0; % keep only brain
        [vT, dT, XT] = sampleT(niiT.hdr, imgT, res, opts.DoF);
        R1 = M * [niiM.hdr.sform_mat; 0 0 0 1];
        hdr = niiM.hdr; hdr.sform_mat = R1(1:3,:);
        msk = nii_xform(opts.mask, hdr);
        F.Values(~logical(msk.img)) = 0; % mask niiM based on current M
        opts.mask = '';
    end
    
    if opts.disp < 3, continue; end
    if ~exist('hsI', 'var') || ~isvalid(hsI(1))
        nii0 = nii_viewer('LocalFunc', 'nii_reorient', niiM, false, 1);
        hs = guidata(nii_viewer(niiT, nii0));
        hsI = findobj(hs.ax(1), 'Type', 'image');
    end
    hsI(1).UserData.Ri = inv(M*[nii0.hdr.sform_mat; 0 0 0 1]);
    nii_viewer('LocalFunc', 'set_cdata', hs); drawnow;
end
if iter>255, warning('nii_coreg:IterExceed', 'Max iterations reached.'); end
mss = std(vT(ind)-vM(ind)*b(1), 1) / muT; % make it independent of img intensity
if opts.disp > 0
    fprintf('iter=%i mss=%4.2f overlap=%.2g\n', iter, mss, mean(ind));
end
if opts.disp == 2
    nii_viewer(niiT, nii_tool('update', niiM, M*[niiM.hdr.sform_mat; 0 0 0 1]));
end
  
%% Return 4x4 transformation matrix from a vector of 3, 6, 7, 9 or 12
function A = affine_mat(P)
% P(1:3)   - xyz translation
% P(4:6)   - xyz rotation about pitch, roll, yaw (degrees)
% P(7)     - single scaling for xyz if numel(P)==7
% P(7:9)   - xyz scaling if numel(P)>=9
% P(10:12) - xyz shearing

A = eye(4); A(1:3, 4) = P(1:3);
if numel(P)==3, return; end
ca = cosd(P(4:6)); sa = sind(P(4:6));
A(1:3, 1:3) = [1 0 0; 0 ca(1) sa(1); 0 -sa(1) ca(1)] * ...
              [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)] * ...
              [ca(3) sa(3) 0; -sa(3) ca(3) 0; 0 0 1];
          
if numel(P)==6, return; end
if numel(P)==7, A(1:3, 1:3) = A(1:3, 1:3) * P(7); return; end
A(1:3, 1:3) = A(1:3, 1:3) * diag(P(7:9));
if numel(P)==9, return; end
S = eye(4); S([5 9 10]) = P(10:12);
A = A * S;

%% Return meaningful voxels, and their derivatives & coordinates
function [vT, dT, XT, RT, p0] = sampleT(hdr, img, res, DoF)
d = hdr.dim(2:4);
dv = res ./ hdr.pixdim(2:4);
[i, j, k] = ndgrid(1:dv(1):d(1)-.5, 1:dv(2):d(2)-.5, 1:dv(3):d(3)-.5);
XT = [i(:)'; j(:)'; k(:)'];
st = rng('default'); XT = XT + rand(size(XT))*0.5; rng(st); % make rst consistant 
XT(4,:) = 1;

RT = [hdr.sform_mat; 0 0 0 1] * [eye(4,3) [-1 -1 -1 1]'];
F = griddedInterpolant(img, 'linear', 'none');
vT = F(XT(1,:), XT(2,:), XT(3,:));
p0 = [0 0 0  0 0 0  1 1 1  0 0 0]'; p0 = p0(1:DoF);
dd = repelem([0.5; 0.5; 1e-4; 1e-2], 3); dd = dd(1:DoF);
dT = zeros(DoF, numel(vT));
XT = RT * XT;
for i = 1:DoF
    p1 = p0; p1(i) = p1(i) + dd(i);
    X = affine_mat(p1) * RT \ XT;
    dT(i,:) = F(X(1,:), X(2,:), X(3,:)) - vT;
end
dT(isnan(dT)) = 0;
dT = dT ./ dd;

a = dT ./ std(dT,1,2); % dT range quite different for shear and zoom
a = var(a);
ind = a > std(a)/50; % arbituray threshold: mainly exclude voxels outside scalp

vT = vT(ind);
dT = dT(:, ind);
XT = XT(:, ind);
%%