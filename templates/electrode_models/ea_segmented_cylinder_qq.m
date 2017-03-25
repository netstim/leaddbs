function [no,fc,seeds,node,elem]=segmented_cylinder(varargin)
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

n1=[r1*sin(t(:)) r1*cos(t(:)) zeros(size(t(:)))];
n2=[r2*sin(t(:)) r2*cos(t(:)) zeros(size(t(:)))];

n3=[r1*sin(t(:)) r1*cos(t(:)) -d-(a*r1*sin(t(:))+b*r1*cos(t(:)))/c];
n4=[r2*sin(t(:)) r2*cos(t(:)) -d-(a*r2*sin(t(:))+b*r2*cos(t(:)))/c];

no=[n1;n2;n3;n4];

%% PLCs of the side-cut fiber

fc={};
count=1;  
for i=1:N-1 % inner and outer cylinder facets
   % the last number in each cell is the fc id
   fc{count}={[i+N i+3*N i+3*N+1 i+N+1],1}; count=count+1;
   fc{count}={[i i+2*N i+2*N+1 i+1],2}; count=count+1;
end
i=N;
fc{count}={[i+N i+3*N 1+3*N 1+N],1}; count=count+1; % inner and outer cylinder facets
fc{count}={[i i+2*N 1+2*N 1],2}; count=count+1; % inner and outer cylinder facets

fc{count}={1:1+N-1,3};count=count+1;  % bottom inner circle
fc{count}={1+N*2:1+N*3-1,5};count=count+1;  % top inner circle

segseed=[];

for i=1:divlen:N-1  % the conductor
    fc{count}={[i i+N i+3*N i+2*N],2}; count=count+1;  % vertical separator 1
    fc{count}={[i+seglen-1 i+seglen-1+N i+seglen-1+3*N i+seglen-1+2*N],2}; count=count+1; % vertical separator 2
    fc{count}={[i:i+seglen-1 i+seglen-1+N:-1:i+N],4};count=count+1; % button outter circle seg
    fc{count}={[i:i+seglen-1 i+seglen-1+N:-1:i+N]+2*N,6};count=count+1;  % top outter circle seg
    segseed=[segseed; [(r1+r2)*[sin(t(i+1)) cos(t(i+1))] h]*0.5];
end
for i=1:divlen:N-divlen % the insulator
    fc{count}={[i+seglen-1:i+divlen i+divlen+N:-1:i+seglen+N-1],7};count=count+1;
    fc{count}={[i+seglen-1:i+divlen i+divlen+N:-1:i+seglen+N-1]+2*N,8};count=count+1;  % top outter circle
end
i=i+divlen;
fc{count}={[i+seglen-1:N 1 1+N 2*N:-1:i+seglen+N-1],7};count=count+1;
fc{count}={[i+seglen-1:N 1 1+N 2*N:-1:i+seglen+N-1]+2*N,8};count=count+1;  % top outter circle

%seeds=[0 0 h*0.5; segseed]; % if you want to label the center cylinder
seeds=segseed;
if(nargout>=4)
    [node,elem,face]=surf2mesh(no,fc,min(no),max(no),1,50,seeds,[],0);
    figure
    plotmesh(node,elem,'x>0 | y>0');
end
