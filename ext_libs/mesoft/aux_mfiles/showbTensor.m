function showbTensor(ten)

%%
for k = 1:size(ten,3), 
    [U D] = eigs(ten(:,:,k)); 
    dir(:,k) = U(:,1)*sqrt(D(1,1)); 
end; 
dir = dir;%.*repmat(sign(dir(3,:)),[3 1]);
figure(1000) ; 
plot3(dir(1,:),dir(2,:),dir(3,:),'*'); hold on;
plot3(-dir(1,:),-dir(2,:),-dir(3,:),'g*'); hold off;