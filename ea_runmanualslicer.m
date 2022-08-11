function [coords_mm,trajectory,markers] = ea_runmanualslicer(options)
%% Function to manually reconstruct electrode using Slicer's fiducial markers
%  Last Revision: 21/05/2018
%  Thushara Perera (c) 2018 Bionics Institute
%  Input:
%   - Lead-DBS options struct
%  Output:
%   - reconstruction data for Lead-DBS pipeline

fprintf(['\nSlicer Instructions\n',...
    '-------------------\n',...
    '1. Move head fiducial marker so it overlays on the most superior contact\n',...
    '2. Move tail fiducial marker so it overlays on the most inferior contact\n',...
    '3. Save the Slicer project overwrite existing files (particularly the fiducial marker file *.fcsv)\n',...
    '4. Exit Slicer\n\n',...
    '--- Alternatively ----\n',...
    '1. Delete all existing markers within Slicer\n',...
    '2. For each side, mark head first, then tail with new fiducials\n',...
    '3. Name of fiducial is not important, but order must always be head first, tail second\n',...
    '4. Save the project overwriting existing fiducial marker (*.fcsv) file\n',...
    '5. Exit Slicer\n',...
    '-------------------\n\n']);
ea_runslicer(options, -1);
disp('returned');
% Read Slicer fiducial markup file
fiducial_path = setBIDSEntity(options.subj.recon.recon, 'desc', 'electrodefiducials', 'ext', 'fcsv');
fid = fopen(fiducial_path, 'r');
F = textscan(fid, '%*s %f %f %f %*[^\n]', 'HeaderLines', 3, 'Delimiter', ',');
fclose(fid);
F = cell2mat(F);

for side = options.sides
    idx = side*2-1;
    markers(side).head = F(idx,:);
    markers(side).tail = F(idx+1,:);

    % add x and y (copied from ea_runmanual.m)
    [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
    markers(side).x = markers(side).head +  xunitv*(options.elspec.lead_diameter/2);
    markers(side).y = markers(side).head + yunitv*(options.elspec.lead_diameter/2);
end

[coords_mm,trajectory,markers] = ea_resolvecoords(markers,options,0);
disp('Manual reconstruction completed.');
