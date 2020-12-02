function ext=ea_getantstransformext(directory)

gz=dir([directory,'glanatComposite.nii.gz']);
h5=dir([directory,'glanatComposite.h5']);

if ~isempty(gz) && isempty(h5)
    ext='.nii.gz';
elseif ~isempty(h5) && isempty(gz)
    ext='.h5';
elseif ~isempty(gz) && ~isempty(h5)
    ext='.nii.gz';
else
    ext='.h5'; % assume default.
end
