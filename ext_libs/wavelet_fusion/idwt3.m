% V = idwt3(W)
%
% Implements the inverse discrete 3D wavelet transform, using Matlab 2D
% functions.
%
% See also: dwt3

% Code written by Iman Aganj.

function V = idwt3(W)

sC = size3(W.dec{1,1,1});
for i=1:3
    s(i) = length(idwt(zeros(sC(i),1), zeros(sC(i),1), W.wl));
end
V = zeros(s);

F = cell(2);
for i=1:2
    for j=1:2
        F{i,j} = zeros([sC(1:2) s(3)]);
    end
end
for i=1:2
    for j=1:2
        for x=1:sC(1)
            for y=1:sC(2)
                F{i,j}(x,y,:) = idwt(W.dec{i,j,1}(x,y,:), W.dec{i,j,2}(x,y,:), W.wl);
            end
        end
    end
end
for k=1:s(3)
    V(:,:,k) = idwt2(F{1,1}(:,:,k), F{2,1}(:,:,k), F{1,2}(:,:,k), F{2,2}(:,:,k), W.wl);
end
