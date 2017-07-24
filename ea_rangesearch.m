function [idx,dist]=rangesearch(c,r,X,mode)
% RANGESEARCH Range Search to find all points within a range.
% [index,distance]=rangesearch(c,r,X,mode) returns all points,x of X which
% are in the range of ||x-c||<r.
% Inputs:
%    c: 1 x d the query vector
%    r: scalar defines the range
%    X: n x d array of all search points
% mode: range mode, 1 - box range, 2 (default) - radial range 
% Outputs:
%    index: indices of points in the range.
% distance: distances between the reference point and points in the range.
%
% See Also: ea_pdist, kdtree, knnsearch

% Version 2.0 by Yi Cao at Cranfield University on 6th April 2008
%

%Examples
%Example 1: Radial range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X);
t=0:0.02:2*pi;
x=c(1)+r*cos(t);
y=c(2)+r*sin(t);
subplot(211)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Radial range search')
%}
%Example 2: Box range 
%{
X=rand(1000,2);
c=[0.5 0.5];
r=0.2;
idx=rangesearch(c,r,X,1);
t=0:0.02:2*pi;
x=c(1)+[-r r r -r -r];
y=c(2)+[r r -r -r r];
subplot(212)
plot(X(:,1),X(:,2),'b.',c(1),c(2),'c+',X(idx,1),X(idx,2),'g*',x,y,'r-','linewidth',2)
title('Box range search')
%}
% Example 3: Large data set search
%{
N=250000;
d=10;
X=randn(N,d);
c=randn(1,d);
r=1;
tic
idx=rangesearch(c,r,X);
toc
%}
% The time is similar to using kdtree (FEX ID 4586).

if nargin<4
    mode=2;
end
if mode==2
    [idx,dist]=range2(c,r,X);
else
    [idx,dist]=range1(c,r,X);
end

function [idx,dist]=range1(c,r,X)
[nPoints,nVariables]=size(X);
s=zeros(nPoints,1);
for d=1:nVariables
    x=abs(X(:,d)-c(d));
    id=x>s;
    s(id)=x(id);
end
fidx=s<r;
idx=find(fidx);
dist=s(fidx);

function [idx,dist]=range2(c,r,X)
nVariables=numel(c);
r2=r*r;
s=0;
for d=1:nVariables
    s=s+(X(:,d)-c(d)).^2;
end
fidx=s<r2;
idx=find(fidx);
dist=s(fidx);

% function [idx,dist]=range1(c,r,X)
% [nPoints,nVariables]=size(X);
% fidx=false(1,nPoints);
% dist=zeros(1,nPoints);
% for k=1:nPoints
%     x=abs(X(k,:)-c);
%     s=0;
%     for d=1:nVariables
%         if x(d)>s
%             s=x(d);
%             if s>r
%                 break
%             end
%         end
%     end
%     if s<r
%         fidx(k)=true;
%         dist(k)=s;
%     end
% end
% idx=find(fidx);
% dist=dist(fidx);

% function [idx,dist]=range2(c,r,X)
% [nPoints,nVariables]=size(X);
% fidx=false(1,nPoints);
% dist=zeros(1,nPoints);
% r2=r*r;
% for k=1:nPoints
%     x=X(k,:)-c;
%     s=0;
%     for d=1:nVariables
%         s=s+x(d)^2;
%     end
%     if s<r2
%         dist(k)=sqrt(s);
%         fidx(k)=true;
%     end
% end
% idx=find(fidx);
% dist=dist(fidx);
