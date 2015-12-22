function ret = fileexists(fname)
ret = false;
fi = fopen(fname);
if fi == -1,
    return;
else
    ret = true;
    fclose(fi);
    return;
end;