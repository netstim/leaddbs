function sinterpstruct = sphereInterpolLUT(opos,sym)


if sym,
    pos = [opos ; -opos];
else
    pos = opos;
end;

N = 150;

beta = 0.5;


s = rescale((1:N)/N,beta);
[x y] = ndgrid([-fliplr(s) 0 s]);
z = sqrt(1-x.^2-y.^2);


X(1,:,:) = x;
Y(1,:,:) = y;
Z(1,:,:) = z;

X(2,:,:) = x;
Y(2,:,:) = y;
Z(2,:,:) = -z;

X(3,:,:) = z;
Y(3,:,:) = y;
Z(3,:,:) = x;

X(4,:,:) = -z;
Y(4,:,:) = y;
Z(4,:,:) = x;

X(5,:,:) = x;
Y(5,:,:) = z;
Z(5,:,:) = y;

X(6,:,:) = x;
Y(6,:,:) = -z;
Z(6,:,:) = y;



pos = pos';


% Compute the convex hull.
C = convhulln(pos');

V = zeros(6,size(X,2),size(X,3));
bary = zeros(3,6,size(X,2),size(X,3));
indices = zeros(3,6,size(X,2),size(X,3));
for i = 1:size(C,1)
    
    %% find points in triangle
   j = C(i,[1 2 3 1]);
   n1 = -cross(pos(:,j(1)) - pos(:,j(2)),pos(:,j(1))); n1 = n1 / norm(n1);
   n2 = -cross(pos(:,j(2)) - pos(:,j(3)),pos(:,j(2))); n2 = n2 / norm(n2);
   n3 = -cross(pos(:,j(3)) - pos(:,j(1)),pos(:,j(3))); n3 = n3 / norm(n3);   
   cos1 = n1(1)*X + n1(2)*Y + n1(3)*Z;
   cos2 = n2(1)*X + n2(2)*Y + n2(3)*Z;
   cos3 = n3(1)*X + n3(2)*Y + n3(3)*Z;
   idx = find(cos1(:)>0 & cos2(:)>0 & cos3(:)>0 & imag(X(:)) == 0 & imag(Y(:)) == 0 & imag(Z(:)) == 0);
   
    %% compute barycentric coordinates with respect to the current triangle
   F = cross(pos(:,j(1))-pos(:,j(2)),pos(:,j(1))-pos(:,j(3))); F = F /norm(F);
   b = F'*pos(:,j(1));
   Pts = [X(idx)' ; Y(idx)' ; Z(idx)'];
   a = F'*Pts;
   Pts = b*Pts./ repmat(a,[3 1]); 
   M = inv(pos(:,j([1 2 3])));   
   bary(:,idx) = M*Pts;
   
    %% set triamngle indices
   indices(:,idx) = repmat(j([1 2 3])',[1 size(idx)]);
   
   %% visualize
%   col = rand;
%   V(idx) = col;      
%   patch(pos(1,j),pos(2,j),pos(3,j),col);

end

if sym,
    pos = opos';
    indices(indices>0) = mod(indices(indices>0) - 1,size(pos,2))+1;
end;




sinterpstruct.barycoords = single(bary);
sinterpstruct.indices = single(indices);
sinterpstruct.beta = single(beta);
sinterpstruct.numpoints = single(size(pos,2));
sinterpstruct.bDir = single(pos);

%testLUTm(sinterpstruct);


% 
% plot3(pos(1,:),pos(2,:),pos(3,:),'ko','markerfacecolor','k');
% 
% % Modify the view.
% vi = [99 36];
% 
% view(vi(1),vi(2)); axis equal
% colormap(spring)
% rotate3d on
% 
% 

%plot3(pos(1,:),pos(2,:),pos(3,:),'ko','markerfacecolor','k');

return;



disp('testing')
n = 200;
[X Y Z] = sphere(n);
C = zeros(n+1,n+1);

for i = 1:(n+1),
    for j = 1:(n+1),
      
            if Z(i,j) > 0.5,
                x = floor((invrescale(X(i,j),beta)+1)*N);
                y = floor((invrescale(Y(i,j),beta)+1)*N);
                C(i,j) = V(1,x,y);
                continue;
            end;
            if Z(i,j) < -0.5,
                x = floor((invrescale(X(i,j),beta)+1)*N);
                y = floor((invrescale(Y(i,j),beta)+1)*N);
                C(i,j) = V(2,x,y);
                continue;
            end;
            if X(i,j) > 0.5,
                z = floor((invrescale(Z(i,j),beta)+1)*N);
                y = floor((invrescale(Y(i,j),beta)+1)*N);
                C(i,j) = V(3,z,y);
                continue;
            end;           
            if X(i,j) < -0.5,
                z = floor((invrescale(Z(i,j),beta)+1)*N);
                y = floor((invrescale(Y(i,j),beta)+1)*N);
                C(i,j) = V(4,z,y);
                continue;
            end;      
            if Y(i,j) > 0,
                x = floor((invrescale(X(i,j),beta)+1)*N);
                z = floor((invrescale(Z(i,j),beta)+1)*N);
                C(i,j) = V(5,x,z);
                continue;
            end;    
            if Y(i,j) <= 0,
                x = floor((invrescale(X(i,j),beta)+1)*N);
                z = floor((invrescale(Z(i,j),beta)+1)*N);
                C(i,j) = V(6,x,z);
                continue;
            end;           
            
        
    end
end;

figure(2);
clf;
surface(X,Y,Z,C,'EdgeColor','none'); hold on;
plot3(pos(1,:),pos(2,:),pos(3,:),'ko','markerfacecolor','k');
view(vi(1),vi(2)); axis equal
colormap(spring)
rotate3d on

% 
% function t = invrescale(f,beta)
% t = f;
% 
% 
% function f = rescale(t,beta)
% f=t;




function t = invrescale(f,beta)
a = 1/(sqrt(1+beta)-sqrt(beta));
b = 1/(1-sqrt(1/beta + 1));
t = sign(f).*(((abs(f)-b)/a).^2-beta);


function f = rescale(t,beta)

a = 1/(sqrt(1+beta)-sqrt(beta));
b = 1/(1-sqrt(1/beta + 1));
f = a*sqrt(beta+t) + b;








