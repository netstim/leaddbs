% Vf = fuseVolW_multi(Vx, Vy, Vz, wl, range, N)
%
% This function uses fuseVolW to run the wavelet-based 3D image fusion
% eight times, each time shifting the volume in one, two, or three of the
% dimensions by (by default) one voxel, and averaging the results. range is
% the possible shifts in each direction (default: [0 1]).
%
% For more information, see also: fuseVolW.
%
% Matlab code written by Iman Aganj. Reference:
%
% I. Aganj, C. Lenglet, E. Yacoub, G. Sapiro, and N. Harel, “A 3D wavelet
% fusion approach for the reconstruction of isotropic-resolution MR images
% from orthogonal anisotropic-resolution scans,” Magnetic Resonance in
% Medicine, vol. 67, no. 4, pp. 1167–1172, 2012.

function Vf = fuseVolW_multi(Vx, Vy, Vz, wl, range, N)

if ~exist('N', 'var')
    N = 1;
end

if ~exist('range', 'var')
    range = [0 1];
end

if ~exist('wl', 'var')
    wl = 'haar';
end

if ~isempty(Vx)
    s = size(Vx);
elseif ~isempty(Vy)
    s = size(Vy);
elseif ~isempty(Vz)
    s = size(Vz);
else
    error('At least one image must be non-empty.')
end

Vf = zeros(s); Vn = Vf; ons = ones(s);
for i = range
    for j = range
        for k = range
            disp([num2str(i) ', ' num2str(j) ', ' num2str(k) ' ...'])
            Vf = Vf + shiftVol(fuseVolW(shiftVol(Vx, i, j, k), shiftVol(Vy, i, j, k), shiftVol(Vz, i, j, k), wl, N), -i, -j, -k);
            Vn = Vn + shiftVol(ons, -i, -j, -k);
        end
    end
end
Vf = Vf ./ Vn;

function V = shiftVol(V0, i, j, k)

if isempty(V0)
    V = [];
else
    V = zeros(size(V0));
    V(1+max(i,0):end+min(i,0), 1+max(j,0):end+min(j,0), 1+max(k,0):end+min(k,0)) = V0(1-min(i,0):end-max(i,0), 1-min(j,0):end-max(j,0), 1-min(k,0):end-max(k,0));
end
