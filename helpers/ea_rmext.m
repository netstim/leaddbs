function fn=ea_rmext(fn)
if iscell(fn)
    [~,fn]=cellfun(@fileparts,fn,'UniformOutput',0); % loose .gz
    [~,fn]=cellfun(@fileparts,fn,'UniformOutput',0); % loose .nii
else
    [~,fn]=fileparts(fn); % loose .gz
    [~,fn]=fileparts(fn); % loose .nii
end