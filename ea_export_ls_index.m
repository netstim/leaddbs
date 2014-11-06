function ea_export_ls_index(serverdir,options)
% very simple static dir2html export function.


fid = fopen([options.prefs.ls.dir,'index.html'], 'w+');

fwrite(fid,'<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">');
fwrite(fid,'<html xmlns="http://www.w3.org/1999/xhtml">');
fwrite(fid,'<head>');
fwrite(fid,'<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />');
fwrite(fid,'<title>LEAD-Server</title>');
fwrite(fid,'<style type="text/css">');
fwrite(fid,'<!--');
fwrite(fid,' body {');
fwrite(fid,'  font-family: Arial;');
fwrite(fid,' }');
fwrite(fid,'//-->');
fwrite(fid,'</style>');
fwrite(fid,'</head>');
fwrite(fid,'');
fwrite(fid,'<body>');
fwrite(fid,'<p><h2>Welcome to Lead Server</h2></p>');

pts=dir([options.prefs.ls.dir,'data']);

for pt=1:length(pts)

    
    if pts(pt).isdir && ~strcmp(pts(pt).name,'.') && ~strcmp(pts(pt).name,'..')
        fwrite(fid,['<p><a href="data/',pts(pt).name,'/index.html">',pts(pt).name,'</a>']);
        
            inpts=dir([options.prefs.ls.dir,'data',filesep,pts(pt).name]);

        for inpt=1:length(inpts)
            
            if inpts(inpt).isdir && ~strcmp(inpts(inpt).name,'.') && ~strcmp(inpts(inpt).name,'..') && ~strcmp(inpts(inpt).name,'slices')
                fwrite(fid,['<br><a href="data/',pts(pt).name,'/',inpts(inpt).name,'/index.html">+    Stimulation  ',inpts(inpt).name,'</a>']);
            elseif inpts(inpt).isdir && strcmp(inpts(inpt).name,'slices')
                fwrite(fid,['<br><a href="data/',pts(pt).name,'/',inpts(inpt).name,'/index.html">+    2D Slice views</a>']);
            end
        end
        
        fwrite(fid,['</p>']);
        
        
        
    end
end
fwrite(fid,'</body>');
fwrite(fid,'</html>');


fclose(fid);

disp('LEAD Server updated.');
    
        disp(sprintf('Enter the following lines to terminal in order to start a test-server: \n'));
        disp(['cd ',serverdir,]);
        disp(['node server.js']);
        disp(sprintf('\n LEAD-Server should then be established successfully.'));
        disp('Navigate to http://localhost:5000/ in your browser to see results.');
    
