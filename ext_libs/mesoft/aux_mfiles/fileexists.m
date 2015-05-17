function ret = fileexists(fname)

if iscell(fname),
    ret = true;
    for k = 1:length(fname),
        ret = ret && fileexists(fname{k});
    end;
    return ;
end;

ret = false;
fi = fopen(fname);
if fi == -1,
    return;
else
    ret = true;
    fclose(fi);
    return;
end;