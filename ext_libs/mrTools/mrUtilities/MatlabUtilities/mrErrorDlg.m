function mrErrorDlg(errstr)
%
% function mrErrorDlg(errstr)
%
% Uses 'verbose' preference to either open a matlab errordlg or just signal
% an error.
%
% To set the 'verbose' preference:
%    mrSetPref('verbose','Yes');
%    mrSetPref('verbose','No');
%
% djh, 5/2005
%
% djh, 5/2007, modified to use mrGetPref instead of Matlab's getpref

verbose = mrGetPref('verbose');

if strcmp(verbose,'Yes')
  errordlg(errstr,'Error!');
  error(errstr);
else
  disp(sprintf('mrERROR: %s',errstr));
  keyboard
end
