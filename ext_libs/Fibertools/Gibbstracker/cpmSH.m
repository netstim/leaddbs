

function SH = cpmSH(dirin,L)

if size(dirin,3) > 1,
    for k = 1:size(dirin,3),
        [V D] = eigs(dirin(:,:,k));
        [~,ix]=sort(D(logical(eye(length(D)))));
        dir(:,k) = V(:,ix(3));

    end;
else
    dir = dirin;
end;

N = size(dir,2);   
costheta = (dir(3,:));
xiy = dir(1,:) + i*dir(2,:);
xiy = xiy ./(abs(xiy)+0.0001);
for k = 1:2:L,
    LegP = legendre(k-1,costheta,'norm');
    for m = -(k-1):(k-1),
        if m < 0 
             SH{k}(m+k,:) = LegP(abs(m)+1,:) .* conj(xiy).^(abs(m)) * (-1)^m;
        else
             SH{k}(m+k,:) = LegP(abs(m)+1,:) .* xiy.^(abs(m));
        end;
    end;  
end;