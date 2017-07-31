% W = dwt3(V, wl)
%
% Implements the discrete 3D wavelet transform, using Matlab 2D functions.
%
% See also: idwt3

% Code written by Iman Aganj.

function W = dwt3(V, wl)

for i=1:3
    s(i) = length(dwt(zeros(size(V,i),1), wl));
end
W.wl = wl;
W.dec = cell(2,2,2);
for i=1:2
    for j=1:2
        for k=1:2
            W.dec{i,j,k} = zeros(s);
        end
    end
end
LL = zeros([s(1:2) size(V,3)]); HL = LL; LH = LL; HH = LL;
for k=1:size(V,3)
    [LL(:,:,k), HL(:,:,k), LH(:,:,k), HH(:,:,k)] = dwt2(V(:,:,k), wl);
end
for x=1:s(1)
    for y=1:s(2)
        [W.dec{1,1,1}(x,y,:), W.dec{1,1,2}(x,y,:)] = dwt(LL(x,y,:), wl);
        [W.dec{2,1,1}(x,y,:), W.dec{2,1,2}(x,y,:)] = dwt(HL(x,y,:), wl);
        [W.dec{1,2,1}(x,y,:), W.dec{1,2,2}(x,y,:)] = dwt(LH(x,y,:), wl);
        [W.dec{2,2,1}(x,y,:), W.dec{2,2,2}(x,y,:)] = dwt(HH(x,y,:), wl);
    end
end
