function err=ea_submitcmd(cmd)

if ~ispc
    err=system(['bash -c "', cmd, '"']);
else
    err=system(cmd);
end
