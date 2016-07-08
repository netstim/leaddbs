 function [Cylinder, EndPlate1, EndPlate2] = ea_cylinder(X1,X2,r,n,cyl_color,closed,lines)
%
% This function constructs a cylinder connecting two center points 
% 
% Usage :
% [Cylinder EndPlate1 EndPlate2] = Cylinder(X1+20,X2,r,n,'r',closed,lines)
%    
% Cylinder-------Handle of the cylinder
% EndPlate1------Handle of the Starting End plate
% EndPlate2------Handle of the Ending End plate
% X1 and X2 are the 3x1 vectors of the two points
% r is the radius of the cylinder
% n is the no. of elements on the cylinder circumference (more--> refined)
% cyl_color is the color definition like 'r','b',[0.52 0.52 0.52]
% closed=1 for closed cylinder or 0 for hollow open cylinder
% lines=1 for displaying the line segments on the cylinder 0 for only
% surface
% 
% Typical Inputs
% X1=[10 10 10];
% X2=[35 20 40];
% r=1;
% n=20;
% cyl_color='b';
% closed=1;
% 


% Calculating the length of the cylinder
length_cyl=norm(abs(X2-X1));

% Creating a circle in the YZ plane
t=linspace(0,2*pi,n)';
x2=r*cos(t);

sintheta = sin(t); sintheta(n) = 0;

x3=r*sintheta;



% Creating the points in the X-Direction
x1=[length_cyl 0];

% Creating (Extruding) the cylinder points in the X-Directions
xx3=repmat(x1,length(x2),1);
xx2=repmat(x2,1,2);
xx1=repmat(x3,1,2);

% Drawing two filled cirlces to close the cylinder
if closed==1
    hold on
    %EndPlate1=fill3(xx1(:,1),xx2(:,1),xx3(:,1),'r');
    
    f.vertices=[xx1(:,1),xx2(:,1),xx3(:,1)];
    f.vertices=[mean(f.vertices);f.vertices];
    f.faces=zeros(length(xx1)-1,3);
    for face=1:length(f.faces)
        f.faces(face,:)=[1,face+1,face+2];
    end
    EndPlate1=patch('Vertices',f.vertices,'Faces',f.faces);
    
    f.vertices=[xx1(:,2),xx2(:,2),xx3(:,2)];
    f.vertices=[[0,0,0];f.vertices];
    f.faces=zeros(length(xx1)-1,3);
    for face=1:length(f.faces)
        f.faces(face,:)=[1,face+1,face+2];
    end
    EndPlate2=patch('Vertices',f.vertices,'Faces',f.faces);
    
    
  
    
end

% Plotting the cylinder along the X-Direction with required length starting
% from Origin
Cylinder=mesh(xx1,xx2,xx3);

% Defining Unit vector along the Z-direction
unit_Vx=[0 0 1];

% Calulating the angle between the x direction and the required direction
% of cylinder through dot product
angle_X1X2=acos( dot( unit_Vx,(X2-X1) )/( norm(unit_Vx)*norm(X2-X1)) )*180/pi;

% Finding the axis of rotation (single rotation) to roate the cylinder in
% X-direction to the required arbitrary direction through cross product
axis_rot=cross([1 0 0],(X2-X1) );

% Rotating the plotted cylinder and the end plate circles to the required
% angles



% if angle_X1X2~=0 % Rotation is not needed if required direction is along X
%     rotate(Cylinder,axis_rot,angle_X1X2,[0 0 0])
%     if closed==1
%         
%         rotate(EndPlate1,axis_rot,angle_X1X2,[0 0 0])
%         rotate(EndPlate2,axis_rot,angle_X1X2,[0 0 0])
%     end
% end

% Till now cylinder has only been aligned with the required direction, but
% position starts from the origin. so it will now be shifted to the right
% position

if closed==1
    vx=get(EndPlate1,'Vertices');
    vx(:,1)=vx(:,1)+X2(1); vx(:,2)=vx(:,2)+X2(2); vx(:,3)=vx(:,3)+X2(3);
    set(EndPlate1,'Vertices',vx);
    vx=get(EndPlate2,'Vertices');
    vx(:,1)=vx(:,1)+X2(1); vx(:,2)=vx(:,2)+X2(2); vx(:,3)=vx(:,3)+X2(3);
    set(EndPlate2,'Vertices',vx);
end
set(Cylinder,'XData',get(Cylinder,'XData')+X2(1))
set(Cylinder,'YData',get(Cylinder,'YData')+X2(2))
set(Cylinder,'ZData',get(Cylinder,'ZData')+X2(3))

% Setting the color to the cylinder and the end plates
set(Cylinder,'FaceColor',cyl_color)
if closed==1
    set([EndPlate1 EndPlate2],'FaceColor',cyl_color,'EdgeColor','none')
else
    EndPlate1=[];
    EndPlate2=[];
end


% If lines are not needed making it disapear
if lines==0
    set(Cylinder,'EdgeAlpha',0)
end