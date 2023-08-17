function ea_clear_xattr
% Use clearXattr to clear xattr of excutables in Lead-DBS and SPM

clearXattr(ea_getearoot, '*maci64');
clearXattr(fileparts(which('spm')), '*maci64');
clearXattr(ea_getearoot, '*.app', 'd');
