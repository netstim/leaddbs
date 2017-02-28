function ea_menuinstall(~,~,cmd)

choice=questdlg(['Please confirm to start downloading dataset: ',cmd],'Download additional data','Proceed','Cancel','Proceed');

if ~strcmp(choice,'Proceed')
    return
end

success=ea_checkinstall(cmd);

if ~(success==-1) % user aborted.
if ~success
    errordlg([cmd,' dataset could not be installed.']);
else
    msgbox([cmd,' successfully installed.']);
end
end
