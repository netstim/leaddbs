function status = ea_submitcmd(cmd)

if ~ispc
    status = system(['bash -c "', cmd, '"']);
else
    status = system(cmd);
end
