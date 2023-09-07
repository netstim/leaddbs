function ea_clear_xattr
% Use clearXattr to clear xattr of excutables in Lead-DBS and SPM

clearXattr(ea_getearoot, '*mac[ia]64');
clearXattr(fileparts(which('spm')), '*mac[ia]64');
clearXattr(ea_getearoot, '*.app', 'd');
