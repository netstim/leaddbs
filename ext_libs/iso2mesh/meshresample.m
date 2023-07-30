function [node,elem]=meshresample(v,f,keepratio)
%
% [node,elem]=meshresample(v,f,keepratio)
%
% resample mesh using CGAL mesh simplification utility
%
% author: Qianqian Fang, <q.fang at neu.edu>
% date: 2007/11/12
%
% input:
%    v: list of nodes
%    f: list of surface elements (each row for each triangle)
%    keepratio: decimation rate, a number less than 1, as the percentage
%               of the elements after the sampling
%
% output:
%    node: the node coordinates of the sampled surface mesh
%    elem: the element list of the sampled surface mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%
% proposed bugfix (for leaddbs use case) by Enrico Opri, 2020.

[node,elem]=domeshsimplify(v,f,keepratio);

if(length(node)==0)
    warning(['Your input mesh contains topological defects, and the ',...
       'mesh resampling utility aborted during processing. Now iso2mesh ',...
       'is trying to repair your mesh with meshcheckrepair. ',...
       'You can also call this manually before passing your mesh to meshresample.'] );
    [vnew,fnew]=meshcheckrepair(v,f);
    [node,elem]=domeshsimplify(vnew,fnew,keepratio);
end
[node,I,J]=unique(node,'rows');
elem=J(elem);
saveoff(node,elem,mwpath('post_remesh.off'));

end

% function to perform the actual resampling
function [node,elem]=domeshsimplify(node,elem,keepratio)

    if true
        %This is just forcing the checkrepair for manifold condition from iso2mesh toolbox (using CGAL tool <cgalsimp2>)

        %Considering that its input is not always manifold in leaddbs (e.g. VATmodel/ea_mesh_electrode)
        %we may have to force the checkrepair first. For now I am assuming
        %<manifold> is not a guaranteed condition before resample (forcing check&repair+meshfix)
        [node,elem]=meshcheckrepair(node,elem);
        [node,elem]=meshcheckrepair(node,elem,'meshfix');

        saveoff(node,elem,mwpath('pre_remesh.off'));
        deletemeshfile(mwpath('post_remesh.off'));
        cmd_str=[' "' ea_getExec(mcpath('cgalsimp2')) '" "' mwpath('pre_remesh.off') '" ' num2str(keepratio) ' "' mwpath('post_remesh.off') '"'];
        fprintf('command: %s\n',cmd_str);%useful for DEBUG
        system(cmd_str);
        [node,elem]=readoff(mwpath('post_remesh.off'));
    else
        %drop in replacement for mesh simplification. Using the mesh
        %simplification available in matlab.
        %The original cgalsimp2 can be slower, and requires the input to be <manifold>

        %Considering that its input is not always manifold in leaddbs (e.g. VATmodel/ea_mesh_electrode)
        %we may have to force the checkrepair first. For now I am assuming
        %<manifold> is not a guaranteed condition before resample (forcing check&repair)
        [node,elem]=meshcheckrepair(node,elem);

        %node is vertex, elem is face
        nfv = reducepatch(node,elem,keepratio);
        elem=nfv.faces;
        node=nfv.vertices;

        % Check/repair to make sure it is manifold (as matlab sometimes
        % reduces it to a non-manifold mesh. Used the matlab reducepatch
        % to reduce the number of dependancies.
        % Consider switching to Surf Ice https://github.com/neurolabusc/surf-ice/
        % ref: http://www.alecjacobson.com/weblog/?p=4444
        % other wrappers at (but non-trivial to compile): https://github.com/alecjacobson/gptoolbox/
        [node,elem]=meshcheckrepair(node,elem);%runs options: dupnode, duplicated, isolated, deep
        %deep should have already been executed by the previous line. In case the toolbox changes, I enforce the removal non-manifold vertices
        [node,elem]=meshcheckrepair(node,elem,'deep');
        [node,elem]=meshcheckrepair(node,elem,'meshfix');
    end
end