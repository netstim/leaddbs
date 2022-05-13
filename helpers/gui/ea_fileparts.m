function [pth,name,ext]=ea_fileparts(pth)
%handling path wrapped with quotes "" (this would happen if ea_path_helper was used on a windows platform (ispc))
%handle also the case of bad termination (e.g. only one quote at the
%beginning, and none at the end)
starts_with_quote=false;
ends_with_quote=false;
if strcmp(pth(1),'"')
    starts_with_quote=true;
    pth(1)=[];
end
if strcmp(pth(end),'"')
    ends_with_quote=true;
    pth(end)=[];
end

[pth, name, ext] = fileparts(pth);

if starts_with_quote || ends_with_quote
    pth=['"' pth '"'];
end
