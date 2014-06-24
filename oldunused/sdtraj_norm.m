function sdtraj_norm
figure
cd
load('trajvectors.mat')

for i=1:size(trajvectors,1)
ntrajvectors(i,:)=trajvectors(i,:)/norm(trajvectors(i,:));
end

%figure, plot3(trajvectors(:,1),trajvectors(:,2),trajvectors(:,3),'r.');

X=ntrajvectors(:,1:3);
X1=X(trajvectors(:,1)<0,:); % right
X2=X(trajvectors(:,1)>0,:); % left
[idx1,ctrs1]=kmeans(X1,1);
[idx2,ctrs2]=kmeans(X2,1);
plot3(X1(:,1),X1(:,2),X1(:,3),'r.','MarkerSize',12)
hold on
plot3(X2(:,1),X2(:,2),X2(:,3),'b.','MarkerSize',12)
plot3(ctrs1(:,1),ctrs1(:,2),ctrs1(:,3),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot3(ctrs1(:,1),ctrs1(:,2),ctrs1(:,3),'ko',...
     'MarkerSize',12,'LineWidth',2)
plot3(ctrs2(:,1),ctrs2(:,2),ctrs2(:,3),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot3(ctrs2(:,1),ctrs2(:,2),ctrs2(:,3),'ko',...
     'MarkerSize',12,'LineWidth',2)
%legend('Left Trajectories','Right Trajectories','Centroids',...
%   'Location','NW')
   


covx1=cov(X1);
for sd=1:4
h1(sd)=plot_gaussian_ellipsoid([ctrs1(:,1),ctrs1(:,2),ctrs1(:,3)],covx1,sd);
set(h1(sd),'facealpha',0.1);
zd=get(h1(sd),'zdata');
set(h1(sd),'cdata',repmat(sd+5,length(zd),length(zd)));
set(h1(sd),'EdgeColor','none');
set(h1(sd),'FaceColor','interp');
set(h1(sd),'FaceLighting','gouraud');

end

a=light('Position',[0,0,100]);
a=light('Position',[0,0,-100]);



covx2=cov(X2);
for sd=1:4
h2(sd)=plot_gaussian_ellipsoid([ctrs2(:,1),ctrs2(:,2),ctrs2(:,3)],covx2,sd);
set(h2(sd),'facealpha',0.1);
zd=get(h2(sd),'zdata');
set(h2(sd),'cdata',repmat(sd,length(zd),length(zd)));
set(h2(sd),'EdgeColor','none');
set(h2(sd),'FaceColor','interp');
set(h2(sd),'FaceLighting','gouraud');
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
%axis([-0.45 0.6, -0.7,0.2]);
   
   