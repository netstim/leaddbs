function slice=ea_minmax(slice)
mn=min(slice(:));

slice=slice-mn;
mx=max(slice(:));
slice=slice/mx;