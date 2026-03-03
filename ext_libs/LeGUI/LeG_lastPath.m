function [path,name,ext] = LeG_lastPath(filterspec,dialogtitle)

savedpath = fullfile(fileparts(which('LeG_lastPath.m')),'LeG_lastPath.txt');

startpath = '';
if exist(savedpath,'file')
    fid = fopen(savedpath,'r'); %read only
    startpath = fscanf(fid,'%s');
    fclose(fid);
end
if ~exist(startpath,'dir') || length(startpath)<=1
    startpath = pwd;
end

[name,path] = uigetfile(filterspec,dialogtitle,startpath,'multiselect','on');
if iscell(name)
    [~,name,ext] = cellfun(@fileparts,fullfile(path,name),'uniformoutput',false);
else
    [~,name,ext] = fileparts(fullfile(path,name));
end

if ~(path==0)
    if ~strcmp(path,'\')
        fid = fopen(savedpath,'w'); %creates file for writing and discards contents
        fprintf(fid,'%s',path);
        fclose(fid);
    end
end

