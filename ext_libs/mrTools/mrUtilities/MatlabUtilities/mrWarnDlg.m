function h = mrWarnDlg(warnstr,verbose)
%
% mrWarnDlg(warnstr,verbose)
%
% Calls Matlab's warndlg or warning depending on verbose preference, or depending on argument verbose ('Yes' or 'No')
%
% To set the 'verbose' preference:
%    mrSetPref('verbose','Yes');
%    mrSetPref('verbose','No');
%
% djh, 5/2005
%
% djh, 5/2007, modified to use mrGetPref instead of Matlab's getpref
%
% $Id$

h = [];
if ieNotDefined('verbose')
  verbose = mrGetPref('verbose');
end

% always display warning on command line
disp(sprintf('Warning: %s',warnstr));
if strcmp(verbose,'Yes')
  h = warndlg(warnstr);
  drawnow;
end
