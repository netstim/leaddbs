function sdtraj
figure
cd
load('trajvectors.mat')
X=trajvectors(:,1:2);
X1=X(X(:,1)<0,:);
X2=X(X(:,1)>0,:);
[idx1,ctrs1]=kmeans(X1,1);
[idx2,ctrs2]=kmeans(X2,1);
plot(X1(:,1),X1(:,2),'r.','MarkerSize',12)
hold on
plot(X2(:,1),X2(:,2),'b.','MarkerSize',12)
plot(ctrs1(:,1),ctrs1(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs1(:,1),ctrs1(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs2(:,1),ctrs2(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs2(:,1),ctrs2(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
%legend('Left Trajectories','Right Trajectories','Centroids',...
%   'Location','NW')
   


covx1=cov(X1);
for sd=1:4
plot_gaussian_ellipsoid([ctrs1(:,1),ctrs1(:,2)],covx1,sd);
end


covx2=cov(X2);
for sd=1:4
plot_gaussian_ellipsoid([ctrs2(:,1),ctrs2(:,2)],covx2,sd);
end

% [stdx1]=std(X1);
%    rectangle('Position',[ctrs1(:,1)-stdx1(1)/2,ctrs1(:,2)-stdx1(2)/2,stdx1(1),stdx1(2)],'Curvature',[1 1])
%    stdx1=stdx1*2;
%    rectangle('Position',[ctrs1(:,1)-stdx1(1)/2,ctrs1(:,2)-stdx1(2)/2,stdx1(1),stdx1(2)],'Curvature',[1 1]) 
%    
%    
% [stdx1]=std(X2);
%    rectangle('Position',[ctrs2(:,1)-stdx1(1)/2,ctrs2(:,2)-stdx1(2)/2,stdx1(1),stdx1(2)],'Curvature',[1 1])
%    stdx1=stdx1*2;
%    rectangle('Position',[ctrs2(:,1)-stdx1(1)/2,ctrs2(:,2)-stdx1(2)/2,stdx1(1),stdx1(2)],'Curvature',[1 1]) 
%    
axis([-0.45 0.6, -0.7,0.2]);
   
   