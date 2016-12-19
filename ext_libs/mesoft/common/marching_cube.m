function [vertVc, triAy, errStr]= marching_cube(dataAy, verbHd)
% function [verteces, triAy, errStr]= marching_cube(dataAy, verbHd)
%   
%
% Bjoern W. Kreher
%
% Medical Physics
% Dept. of  Diagnostic Radiology 
% University Hospital Freiburg
%
% 08/04
%
% UNIX


persistent marchC

if isempty(marchC)
    marchC= open('marching_cube.mat');
end

vertVc= []; triAy= []; errStr= '';

if ~exist('verbHd') || ~ishandle(verbHd)
    verbHd= [];
end




mask= zeros(size(dataAy));
mask(dataAy ~= 0)= 1;

dataIndex= 1 + ...
    001*mask(1:(end - 1), 1:(end - 1), 1:(end - 1)) + 002*mask(2:end,       1:(end - 1), 1:(end - 1)) + ...
    004*mask(2:end,       2:end,       1:(end - 1)) + 008*mask(1:(end - 1), 2:end,       1:(end - 1)) + ...
    016*mask(1:(end - 1), 1:(end - 1), 2:end      ) + 032*mask(2:end,       1:(end - 1), 2:end)       + ...
    064*mask(2:end,       2:end,       2:end      ) + 128*mask(1:(end - 1), 2:end,       2:end);


sizeAy= size(dataIndex);
idx= find((dataIndex > 1) & (dataIndex < 256));
% Reserviere Speicher
triAy= zeros(sum(marchC.triNo(dataIndex(idx))), 3);
vertVc= zeros(sum(marchC.vertNo(dataIndex(idx))), 3);


posVc= reshape_index(idx, sizeAy) + 0.5*ones(length(idx), 3);
curVert= 0;     curTri= 0;
for i= 1:length(idx)
    triIdx= marchC.triTable(dataIndex(idx(i)), :);
    k= 1;
    tmp= unique(triIdx);
    vertIdx= tmp(2:end);
    mapIndex= zeros(2, 12);
    mapIndex(vertIdx)= 1:length(vertIdx);
    vertVc((curVert + 1):(curVert + length(vertIdx)), :)= ones(length(vertIdx), 1)*posVc(i, :) + marchC.edgeMidTable(vertIdx, :);
    while (triIdx(k) ~= -1)
        curTri= curTri + 1;
        triAy(curTri, :)= mapIndex(triIdx(k:(k + 2))) + curVert;
        k= k + 3;
    end    
    curVert= curVert + length(vertIdx);
    if mod(i, 100) == 0
        verbStr= sprintf('Reconstruct 3D structure (%g%%)', round(1000*i/length(idx))/10);
        if ~isempty(verbHd)
            set(verbHd, 'String', verbHd);
            drawnow;
        else
            disp(verbStr);
        end
    end
end


%figure, triHd= trisurf(triAy, vertVc(:, 1), vertVc(:, 2), vertVc(:, 3)); axis equal;
%set(triHd, 'LineStyle', 'none', 'FaceLighting', 'phong', 'EdgeLighting', 'phong');
%set(triHd, 'LineStyle', 'none', 'FaceLighting', 'flat', 'EdgeLighting', 'flat');
