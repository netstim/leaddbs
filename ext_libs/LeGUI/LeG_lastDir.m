function seldir = LeG_lastDir(dialogtitle)

saveddir = fullfile(fileparts(which('LeG_lastDir.m')),'LeG_lastDir.txt');

startdir = '';
if exist(saveddir,'file')
    fid = fopen(saveddir,'r'); %read only
    startdir = fscanf(fid,'%s');
    fclose(fid);
end
if ~exist(startdir,'dir') || length(startdir)<=1
    startdir = pwd;
end

seldir = uigetdir(startdir,dialogtitle);
if ~(seldir==0)
    fid = fopen(saveddir,'w'); %creates file for writing and discards contents
    fprintf(fid,'%s',seldir);
    fclose(fid);
end