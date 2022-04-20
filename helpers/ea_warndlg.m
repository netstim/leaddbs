function ea_warndlg(varargin)
% Warper of warndlg to show warning in both dialog and command window.

warning('off', 'backtrace');
warning(sprintf(varargin{:}));
warndlg(sprintf(varargin{:}), 'Warning');
warning('on', 'backtrace');
