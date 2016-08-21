function ext=ea_getantstransformext(directory,options)

gz=dir([directory,'gl',options.prefs.prenii_unormalized,'Composite.nii.gz']);
h5=dir([directory,'gl',options.prefs.prenii_unormalized,'Composite.h5']);

if ~isempty(gz) && isempty(h5)
    ext='.gz';
elseif ~isempty(h5) && isempty(gz)
    ext='.h5';
elseif ~isempty(gz) && ~isempty(h5)
    ext='.gz';
else
    ext='.h5'; % assume default.
end
