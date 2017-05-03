function [node,elem,face]=s2m(v,f,keepratio,maxvol,method,regions,holes)
%
% [node,elem,face]=s2m(v,f,keepratio,maxvol,method)
% [node,elem,face]=s2m(v,f,keepratio,maxvol,'tetgen',regions,holes)
%
% volumetric mesh generation from a closed surface, shortcut for surf2mesh
%
% author: Qianqian Fang (fangq <at> nmr.mgh.harvard.edu)
%
% inputs and outputs are similar to those defined in surf2mesh
%
% if method='cgalpoly', s2m will call cgals2m and keepratio should be a 
% structure (as the 'opt' input in cgals2m)
%
% input default values:
%       method: if ignored, iso2mesh uses surf2mesh ('tetgen') to do the
%               tetrahedral mesh generation
%       regions,holes: if ignored, iso2mesh assumes both are empty
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

if(nargin>=5 && strcmp(method,'cgalpoly'))
    [node,elem,face]=cgals2m(v,f,keepratio,maxvol);
    return;
end
if(nargin<=5)
    regions=[];
end
if(nargin<=6)
    holes=[];
end

[node,elem,face,success]=surf2mesh(v,f,min(v,[],1),max(v,[],1),keepratio,maxvol,regions,holes);

if ~success % try with gmsh next
    error('Something went wrong / probably self-intersecting faces.');
    %     verts=v;
    %     faces=zeros(length(f),4);
    facc=cell2mat(f);
    fv.vertices=v;
    fv.faces=f;
    faaa=patch('faces',f,'vertices',v);
    % convert to .gmsh
    gmsh_mesh3d_write(mwpath('post_vmesh.geo'),3,size(verts,1),verts',4,size(faces,1),faces');
end
