% Vf = fuseVolW(Vx, Vy, Vz, wl, N)
%
% This function implements a wavelet-based 3D image fusion. It takes as
% input three volumes (of the same size) Vx, Vy, and Vz, which have low
% resolution in the first, second, and third dimensions, respectively, and
% high resolution in the two other dimensions. wl is the wavelet basis
% (default: 'haar'), and N is the decomposition level (default: 1). The
% output of the algorithm, Vf, is the fused volume. This code requires the
% Matlab Wavelet Toolbox.
%
% See also: fuseVolW_multi
%
% Matlab code written by Iman Aganj. Reference:
%
% I. Aganj, C. Lenglet, E. Yacoub, G. Sapiro, and N. Harel, “A 3D wavelet
% fusion approach for the reconstruction of isotropic-resolution MR images
% from orthogonal anisotropic-resolution scans,” Magnetic Resonance in
% Medicine, vol. 67, no. 4, pp. 1167–1172, 2012.

function Vf = fuseVolW(Vx, Vy, Vz, wl, N)

V{1} = Vx; V{2} = Vy; V{3} = Vz; clear Vx Vy Vz

if ~exist('wl', 'var')
    wl = 'haar';
end
C0 = cell(3,1);
for i=1:3
    if isempty(V{i})
        C0{i} = cell(2,2,2);
    else
        disp(['Computing the wavelet transform of volume ' num2str(i) ' ...'])
        W = dwt3(V{i}, wl);
        C0{i} = W.dec;
        s0 = size3(V{i});
        s = size3(C0{i}{1,1,1});
    end
end

if ~exist('s', 'var')
    Vf = [];
    return
end

C = cell(2,2,2);
for x = 1:2
    for y = 1:2
        for z = 1:2
            if exist('N', 'var') && N>1
                %if x+y+z==5, C{x,y,z} = C0{12-x-y*2-z*3}{x,y,z};, else
                if x==2
                    C0{1}{x,y,z} = [];
                end
                if y==2
                    C0{2}{x,y,z} = [];
                end
                if z==2
                    C0{3}{x,y,z} = [];
                end
                C{x,y,z} = fuseVolW(C0{1}{x,y,z}, C0{2}{x,y,z}, C0{3}{x,y,z}, wl, N-1);
                if isempty(C{x,y,z})
                    C{x,y,z} = zeros(s);
                end
            else
                C{x,y,z} = zeros(s);
                k = 0;
                for i = 1:3
                    if ~isempty(V{i}) && ((i==1 && x==1) || (i==2 && y==1) || (i==3 && z==1))
                        k = k + 1;
                        C{x,y,z} = C{x,y,z} + C0{i}{x,y,z};
                    end
                end
                if k > 0
                    C{x,y,z} = C{x,y,z} / k;
                end
                %C{2,2,2} = (C0{1}{2,2,2}+C0{2}{2,2,2}+C0{3}{2,2,2})/3;
            end
        end
    end
end

disp('Computing the inverse wavelet transform...')
W.dec = C;
Vf = idwt3(W);
Vf = Vf(1:s0(1), 1:s0(2), 1:s0(3));
