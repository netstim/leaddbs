function ea_exportfiducials(options, filename)
% Save fiducials as comma separated value file in patient folder
%  Last Revision: 10/04/2018
%  Thushara Perera (c) 2018 Bionics Institute
%  Input:
%   - lead dbs options struct
%   - filename with extension *.fcsv for Slicer or *.csv for Excel
%  Output:
%   - fiducial marker file will be saved in reconstruction folder

[coords,~,~] = ea_load_reconstruction(options);
header = ['# Markups fiducial file version = 4.7\r\n',...
          '# CoordinateSystem = 0\r\n',...
          '# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\r\n'];
c = vertcat(coords{:});
fid = fopen(filename, 'w');
fprintf(fid, header);
for i = 1:length(c)
    idx = num2str(i-1);
    fprintf(fid, ['vtkMRMLMarkupsFiducialNode_', idx, ',',num2str(c(i,1)),',',num2str(c(i,2)),',',num2str(c(i,3)),...
        ',0,0,0,1,1,1,0,E', idx, ',,vtkMRMLScalarVolumeNode2\r\n']);
end
fclose(fid);
