function ea_peakcolor(n3)
% shows RGB values in nx3 array

X=zeros(1,size(n3,1),3);
for entry=1:size(n3,1);
    for rgb=1:3
        X(1,entry,rgb)=n3(entry,rgb);
    end
end

figure, imagesc(X);
