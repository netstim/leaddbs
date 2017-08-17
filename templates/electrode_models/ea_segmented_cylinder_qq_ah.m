function [no,fc,seeds,fcc,noc,fci,noi]=segmented_cylinder(varargin)
% [no,fc,seeds]=segmented_cylinder();
%   or
% [no,fc,seeds,node,elem]=segmented_cylinder(ndiv,h,r1,r2,nseg,fillratio);
%
% input:
%    ndiv: angular division (default: 60)
%    h: height of the cylinder (default: 100)
%    r1,r2: inner and outer radii of the cylinder (default: 20,25)
%    nseg: electrode segment number (default 3)
%    fillratio: a ratio between 0-1 (default 0.8)
%
% author: qianqian fang <q.fang at neu.edu>
%

%% angular resolution
ndiv=60;
dt=pi/ndiv;
t=0:dt:2*pi-dt;
%t=linspace(0,2*pi,120);
%t=0:dt:2*pi;

%t2=(0:ndiv)/ndiv*2*pi;

%% parameters of the electrode
h=1;  % length
r1=0.5;  % inner radius
r2=1;  % outer radius
nseg=3; % how many segments
fillratio=0.8; % elctrode filling ratio per segment, 1-solid cylinder

if(nargin>=1)
    ndiv=varargin{1};
elseif(nargin>=2)
    h=varargin{2};
elseif(nargin>=3)
    r1=varargin{3};
elseif(nargin>=4)
    r2=varargin{4};
elseif(nargin>=5)
    nreg=varargin{5};
elseif(nargin>=6)
    fillratio=varargin{6};
end

N=length(t);
a=0;
b=0;
c=1;
d=-h;

divlen=2*pi/nseg/dt;          % nodes per division=conductor+insulator
seglen=2*pi/nseg*fillratio/dt;% nodes per conductor
%% key nodes of a side-cut fiber

n1=[r1*cos(t(:)) r1*sin(t(:)) zeros(size(t(:)))];
n2=[r2*cos(t(:)) r2*sin(t(:)) zeros(size(t(:)))];

n3=[r1*cos(t(:)) r1*sin(t(:)) -d-(a*r1*sin(t(:))+b*r1*cos(t(:)))/c];
n4=[r2*cos(t(:)) r2*sin(t(:)) -d-(a*r2*sin(t(:))+b*r2*cos(t(:)))/c];

no=[n1;n2;n3;n4];

%% PLCs of the side-cut fiber

fc={};
fcc={};
fci={};
count=1;  icount=1; ccount=1;
for i=1:N-1 % inner and outer cylinder facets
    % the last number in each cell is the fc id
    fc{count}=[i+N i+3*N i+3*N+1 i+N+1]; count=count+1;
    fc{count}=[i i+2*N i+2*N+1 i+1]; count=count+1;
end
for j=1:divlen:N % the insulator
    for d=1:(divlen-seglen+1)
        i=j+d+seglen-2;
        fci{icount}=wrapnodes([i+N i+3*N i+3*N+1 i+N+1],N);icount=icount+1;
    end
end


for j=1:divlen:N-1 % the contact
    for d=1:seglen-1
        i=j+d-1;
        fcc{ccount}=[i+N i+3*N i+3*N+1 i+N+1]; ccount=ccount+1;
        fcc{ccount}=[i i+2*N i+2*N+1 i+1]; ccount=ccount+1;
        % add these to insulator, too:
    end
end

for j=1:divlen:N-1 % the insulator
    for d=1:seglen-1
        i=j+d-1;
        % add these to insulator, too:
        fci{icount}=[i i+2*N i+2*N+1 i+1]; icount=icount+1;

    end
end


i=N;
fc{count}=[i+N i+3*N 1+3*N 1+N]; count=count+1; % inner and outer cylinder facets
fc{count}=[i i+2*N 1+2*N 1]; count=count+1; % inner and outer cylinder facets

%fci{icount}=[i+N i+3*N 1+3*N 1+N]; icount=icount+1; % inner and outer cylinder facets
%fci{icount}=[i i+2*N 1+2*N 1]; icount=icount+1; % inner and outer cylinder facets


fc{count}=1:1+N-1;count=count+1;  % bottom inner circle

fc{count}=1+N*2:1+N*3-1;count=count+1;  % top inner circle

segseed=[];

for i=1:divlen:N-1  % the conductor
    fc{count}=[i i+N i+3*N i+2*N]; count=count+1;  % vertical separator 1
    fc{count}=[i+seglen-1 i+seglen-1+N i+seglen-1+3*N i+seglen-1+2*N]; count=count+1; % vertical separator 2
    fc{count}=[i:i+seglen-1 i+seglen-1+N:-1:i+N];count=count+1; % button outter circle seg
    fc{count}=[i:i+seglen-1 i+seglen-1+N:-1:i+N]+2*N;count=count+1;  % top outter circle seg
    segseed=[segseed; [(r1+r2)*[sin(t(i+1)) cos(t(i+1))] h]*0.5];
    
    
    fcc{ccount}=[i i+N i+3*N i+2*N]; ccount=ccount+1;  % vertical separator 1
    fcc{ccount}=[i+seglen-1 i+seglen-1+N i+seglen-1+3*N i+seglen-1+2*N]; ccount=ccount+1; % vertical separator 2
    fcc{ccount}=[i:i+seglen-1 i+seglen-1+N:-1:i+N];ccount=ccount+1; % button outter circle seg
    fcc{ccount}=[i:i+seglen-1 i+seglen-1+N:-1:i+N]+2*N;ccount=ccount+1;  % top outter circle seg
    
    fci{icount}=[i i+N i+3*N i+2*N]; icount=icount+1;  % vertical separator 1
    fci{icount}=[i+seglen-1 i+seglen-1+N i+seglen-1+3*N i+seglen-1+2*N]; icount=icount+1; % vertical separator 2
end
for i=1:divlen:N-divlen % the insulator
    fc{count}=[i+seglen-1:i+divlen i+divlen+N:-1:i+seglen+N-1];count=count+1;
    fc{count}=[i+seglen-1:i+divlen i+divlen+N:-1:i+seglen+N-1]+2*N;count=count+1;  % top outter circle
    
end
i=i+divlen;
fc{count}=[i+seglen-1:N 1 1+N 2*N:-1:i+seglen+N-1];count=count+1;
fc{count}=[i+seglen-1:N 1 1+N 2*N:-1:i+seglen+N-1]+2*N;count=count+1;  % top outter circle



%figure, plotmesh(no,fci)
%keyboard

%figure, plot3(no(:,1),no(:,2),no(:,3),'r.')

%seeds=[0 0 h*0.5; segseed]; % if you want to label the center cylinder
seeds=segseed;
% if(nargout>=4)
%     [node,elem,face]=surf2mesh(no,fc,min(no),max(no),1,50,seeds,[],0);
%     figure
%     plotmesh(node,elem,'x>0 | y>0');
% end
% figure, plotmesh(no,fci)



% custom top/endplates for insulator

% kreisdeckel
% fci{icount}=1:1+N-1;icount=icount+1;  % bottom inner circle  <----
% fci{icount}=1+N*2:1+N*3-1;icount=icount+1;  % top inner circle <----
% 
% for i=1:divlen:N-divlen % the insulator
%     
%     fci{icount}=[i+seglen-1:i+divlen i+divlen+N:-1:i+seglen+N-1];icount=icount+1; % <--
%     fci{icount}=[i+seglen-1:i+divlen i+divlen+N:-1:i+seglen+N-1]+2*N;icount=icount+1;  % top outter circle <---
% end

%fci{icount}=[i+seglen-1:N 1 1+N 2*N:-1:i+seglen+N-1];icount=icount+1; % <---
% fci{icount}=[i+seglen-1:N 1 1+N 2*N:-1:i+seglen+N-1]+2*N;icount=icount+1; % <--- top outter circle 

spacelen=(divlen-seglen);

seglenp=seglen+1;
top=[1:seglenp-1,... % inner 1st segment
       (seglenp-1:seglenp+spacelen)+N,... % outer 1st segment
    (seglenp+spacelen):(2*seglenp+spacelen-2),... % inner 2nd segment
    (2*seglenp+spacelen-2)+N:(2*seglenp+2*spacelen)+N-1,... % outer 2nd segment
    (2*seglenp+2*spacelen)-1:(3*seglenp+2*spacelen)-3,... % inner 3rd segment
    (3*seglenp+2*spacelen)+N-3:(3*seglenp+3*spacelen)+N-2]; % outer 3rd segment

top(top>2*N)=top(top>2*N)-N;

bottom=top+2*N;

fci=[fci,{top},{bottom}];

keyboard
% figure, plotmesh(no,fci,'EdgeColor','g','FaceColor','none')
% hold on
% plotmesh(no,fcc,'EdgeColor','r','FaceColor','none')
% 
% plotmesh(no,fc,'EdgeColor','b','FaceColor','none')


%[noi,fci]=removedupnodes(no,fci,1e-6);


%

%[noi,fci]=removeisolatednode(no,fci);

[noi,velem,fci]=s2m(noi,fci,1,3);


[noc,fcc]=removeisolatednode(no,fcc);

[noc,velem,fcc]=s2m(noc,fcc,1,3);

%figure, plotmesh(no,fcc);
%figure,plotmesh(no,fci);
%hold on
%plotmesh(no,fc,'edgecolor','r','facecolor','none');
            
            

function fc=wrapnodes(elem, count)

fc=elem;
idx=find(elem>[2 4 4 2]*count);
fc(idx)=fc(idx)-count;


