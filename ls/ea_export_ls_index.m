function ea_export_ls_index(options)
% very simple static dir2html export function.


fid = fopen([options.prefs.ls.dir,'index.html'], 'wt');

fwrite(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
fwrite(fid,'<html xmlns="http://www.w3.org/1999/xhtml">');
fwrite(fid,'<head>');
fwrite(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />');
fwrite(fid,'<title>Untitled Document</title>');
fwrite(fid,'</head>');
fwrite(fid,'');
fwrite(fid,'<body>');
fwrite(fid,'<p>Welcome to Lead Server</p>');

pts=dir([options.prefs.ls.dir,'data']);
for pt=1:length(pts)
    if pts(pt).isdir && ~strcmp(pts(pt).name,'.') && ~strcmp(pts(pt).name,'..')
        fwrite(fid,['<p><a href="data/',pts(pt).name,'/index.html">',pts(pt).name,'</a></p>']);
    end
end
fwrite(fid,'</body>');
fwrite(fid,'</html>');


fclose(fid);
