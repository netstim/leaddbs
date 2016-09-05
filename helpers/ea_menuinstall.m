function ea_menuinstall(~,~,cmd,force)

choice=questdlg(['Please confirm to start downloading dataset: ',cmd],'Download additional data','Proceed','Cancel','Proceed');

if ~strcmp(choice,'Proceed')
    return
end

success=ea_checkinstall(cmd,force);
if ~success
    errordlg([cmd,' dataset could not be installed.']);
else
    msgbox([cmd,' successfully installed.']);
end